# SPDX-FileCopyrightText:  gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

"""
Compose the Great Britain focused PyPSA network described in
``doc/gb-model/index.rst``.

The rule assembles the clustered PyPSA-Eur base network with GB-specific
artefacts (manual region shapes, neighbouring countries, adjusted grid
connection costs) so that downstream rules can import a consistent
``networks/composed_{clusters}.nc`` snapshot.
"""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import pandas as pd
import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.add_electricity import (
    attach_hydro,
    attach_wind_and_solar,
    load_and_aggregate_powerplants,
    load_costs,
    sanitize_carriers,
    sanitize_locations,
)
from scripts.prepare_sector_network import add_electricity_grid_connection

logger = logging.getLogger(__name__)


@dataclass
class CompositionContext:
    resources_root: Path
    countries: tuple[str, ...]
    costs_path: Path | None
    costs_config: dict[str, Any]
    max_hours: dict[str, Any] | None


def create_context(network_path: str, config: dict[str, Any]) -> CompositionContext:
    resources_root = Path(network_path).parents[1]
    costs_cfg = copy.deepcopy(config.get("costs", {}))
    cost_year = costs_cfg.get("year")
    costs_path = (
        resources_root / f"costs_{cost_year}.csv" if cost_year is not None else None
    )

    return CompositionContext(
        resources_root=resources_root,
        countries=tuple(config.get("countries", [])),
        costs_path=costs_path,
        costs_config=costs_cfg,
        max_hours=config.get("electricity", {}).get("max_hours"),
    )


def _load_powerplants(
    powerplants_path: Path | None,
    costs,
    config: dict[str, Any],
) -> pd.DataFrame:
    clustering_cfg = config.get("clustering", {})
    consider_efficiency = clustering_cfg.get("consider_efficiency_classes", False)
    aggregation_strategies = clustering_cfg.get("aggregation_strategies", {})
    exclude_carriers = clustering_cfg.get("exclude_carriers", [])

    columns = ["bus", "carrier", "p_nom", "max_hours"]

    if powerplants_path is None:
        logger.info("No power plant input provided; continuing without existing capacities")
        return pd.DataFrame(columns=columns)

    if not powerplants_path.exists():
        logger.warning(
            "Power plant data %s not found; continuing without existing capacities",
            powerplants_path,
        )
        return pd.DataFrame(columns=columns)

    return load_and_aggregate_powerplants(
        str(powerplants_path),
        costs,
        consider_efficiency_classes=consider_efficiency,
        aggregation_strategies=aggregation_strategies,
        exclude_carriers=exclude_carriers,
    )


def integrate_renewables(
    network: pypsa.Network,
    config: dict[str, Any],
    costs,
    inputs: Mapping[str, str],
    powerplants_path: Path | None,
    hydro_capacities_path: Path | None,
) -> None:
    electricity_cfg = config.get("electricity", {})
    renewable_carriers = list(electricity_cfg.get("renewable_carriers", []))

    if not renewable_carriers:
        logger.info("No renewable carriers configured; skipping integration")
        return

    if costs is None:
        logger.warning("Cost data unavailable; skipping renewable integration")
        return

    extendable_carriers = copy.deepcopy(electricity_cfg.get("extendable_carriers", {}))
    extendable_carriers.setdefault("Generator", [])

    ppl = _load_powerplants(powerplants_path, costs, config)

    profile_paths = {
        key: Path(value)
        for key, value in inputs.items()
        if key.startswith("profile_")
    }

    available_non_hydro = [
        carrier
        for carrier in renewable_carriers
        if carrier != "hydro" and f"profile_{carrier}" in profile_paths
    ]

    missing_profiles = sorted(
        carrier
        for carrier in renewable_carriers
        if carrier != "hydro" and f"profile_{carrier}" not in profile_paths
    )
    if missing_profiles:
        logger.warning(
            "Renewable profiles missing for carriers %s; skipping them",
            ", ".join(missing_profiles),
        )

    if available_non_hydro:
        landfall_lengths = {
            tech: settings.get("landfall_length")
            for tech, settings in config.get("renewable", {}).items()
            if isinstance(settings, dict) and settings.get("landfall_length") is not None
        }
        line_length_factor = config.get("lines", {}).get("length_factor", 1.0)

        attach_wind_and_solar(
            network,
            costs,
            ppl,
            {k: str(v) for k, v in profile_paths.items()},
            available_non_hydro,
            extendable_carriers,
            line_length_factor,
            landfall_lengths,
        )

    if "hydro" not in renewable_carriers:
        return

    hydro_profile = profile_paths.get("profile_hydro")
    if hydro_profile is None or not hydro_profile.exists():
        logger.warning("Hydro profile not available; skipping hydro integration")
        return

    if hydro_capacities_path is None or not hydro_capacities_path.exists():
        logger.warning(
            "Hydro capacities file %s missing; skipping hydro integration",
            hydro_capacities_path,
        )
        return

    hydro_cfg = copy.deepcopy(config.get("renewable", {}).get("hydro", {}))
    carriers = hydro_cfg.pop("carriers", [])

    attach_hydro(
        network,
        costs,
        ppl,
        str(hydro_profile),
        str(hydro_capacities_path),
        carriers,
        **hydro_cfg,
    )

def add_gb_components(
    n: pypsa.Network,
    context: CompositionContext,
) -> pypsa.Network:
    if context.countries:
        keep = n.buses.country.isin(context.countries)
        drop = n.buses.index[~keep]
        if len(drop) > 0:
            logger.info("Removing %d buses outside target countries", len(drop))
            n.mremove("Bus", drop)

    meta = n.meta.setdefault("gb_model", {})
    if context.countries:
        meta["countries"] = list(context.countries)

    return n


def add_pypsaeur_components(
    n: pypsa.Network,
    config: dict[str, Any],
    context: CompositionContext,
    costs,
) -> pypsa.Network:
    if costs is not None:
        add_electricity_grid_connection(n, costs)
        n.meta.setdefault("gb_model", {})["costs_path"] = str(context.costs_path)

    sanitize_locations(n)
    try:
        sanitize_carriers(n, config)
    except KeyError as exc:  # pragma: no cover - tolerate partial configs
        logger.debug("Skipping carrier sanitisation due to missing config: %s", exc)
    return n


def finalise_composed_network(
    n: pypsa.Network,
    context: CompositionContext,
) -> pypsa.Network:
    n.determine_network_topology()
    meta = n.meta.setdefault("gb_model", {})
    meta["resources_root"] = str(context.resources_root)
    meta["composed"] = True
    n.consistency_check()
    return n


def compose_network(
    network_path: str,
    output_path: str,
    config: dict[str, Any],
    inputs: Mapping[str, str],
) -> None:
    network = pypsa.Network(network_path)
    context = create_context(network_path, config)
    if "tech_costs" in inputs:
        context.costs_path = Path(inputs["tech_costs"])
    add_gb_components(network, context)

    costs = None
    costs_path = context.costs_path
    if costs_path is not None and costs_path.exists():
        weights = getattr(network.snapshot_weightings, "objective", None)
        nyears = float(weights.sum()) / 8760.0 if weights is not None else 1.0
        try:
            costs = load_costs(
                str(costs_path),
                context.costs_config,
                max_hours=context.max_hours,
                nyears=nyears,
            )
        except Exception as exc:  # pragma: no cover - keep composing resilient
            logger.warning(
                "Failed to load cost data from %s: %s",
                costs_path,
                exc,
            )
    elif costs_path is not None:
        logger.info(
            "Cost data %s not found; skipping cost adjustments",
            costs_path,
        )

    powerplants_path = Path(inputs["powerplants"]) if "powerplants" in inputs else None
    hydro_capacities_path = (
        Path(inputs["hydro_capacities"]) if "hydro_capacities" in inputs else None
    )

    integrate_renewables(
        network,
        config,
        costs,
        inputs,
        powerplants_path,
        hydro_capacities_path,
    )

    add_pypsaeur_components(network, config, context, costs)
    finalise_composed_network(network, context)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    network.export_to_netcdf(output_path)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("compose_network", clusters=100)

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    compose_network(
        network_path=snakemake.input.network,
        output_path=snakemake.output.network,
        config=snakemake.config,
        inputs={k: str(v) for k, v in snakemake.input.items()},
    )
