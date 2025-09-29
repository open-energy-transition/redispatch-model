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
from typing import Any

import pypsa

from scripts._helpers import configure_logging, set_scenario_config
from scripts.add_electricity import (
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
    network_path: str, output_path: str, config: dict[str, Any]
) -> None:
    network = pypsa.Network(network_path)
    context = create_context(network_path, config)
    add_gb_components(network, context)

    costs = None
    if context.costs_path is not None and context.costs_path.exists():
        weights = getattr(network.snapshot_weightings, "objective", None)
        nyears = float(weights.sum()) / 8760.0 if weights is not None else 1.0
        try:
            costs = load_costs(
                str(context.costs_path),
                context.costs_config,
                max_hours=context.max_hours,
                nyears=nyears,
            )
        except Exception as exc:  # pragma: no cover - keep composing resilient
            logger.warning(
                "Failed to load cost data from %s: %s",
                context.costs_path,
                exc,
            )
    elif context.costs_path is not None:
        logger.info(
            "Cost data %s not found; skipping cost adjustments",
            context.costs_path,
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
    )
