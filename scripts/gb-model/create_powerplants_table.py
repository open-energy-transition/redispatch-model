# SPDX-FileCopyrightText:  gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


"""
Capacity table generator.

This is a script to GB/Eur capacities defined for / by the FES to fix `p_nom` in PyPSA-Eur.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def _map_names(
    series: pd.Series, mapping: dict[str, str | dict], search_columns: list
) -> str | None:
    """Map carriers/sets to a standard set."""
    try:
        mapped: str | dict | None = mapping.get(series[search_columns[0]])
    except IndexError:
        breakpoint()
    if isinstance(mapped, dict):
        mapped = _map_names(series, mapped, search_columns[1:])
    return mapped


def table_gb_capacities(
    df: pd.DataFrame,
    carrier_mapping: dict,
    set_mapping: dict,
) -> pd.DataFrame:
    """
    To table the powerplant capacities in a format required by PyPSA-Eur

    Args:
        df (pd.DataFrame): FES GSP-level powerplant data table
        carrier_mapping (dict): dictionary to map technologies to PyPSA-Eur carriers names
        set_mapping (dict): dictionary to map technologies to powerplantmatching set names
    """
    search_columns = ["Technology", "Technology Detail"]
    df_cleaned = df.where(df.data > 0).dropna(subset=["data", "Latitude", "Longitude"])
    df_cleaned["carrier"] = df_cleaned.apply(
        _map_names, axis=1, mapping=carrier_mapping, search_columns=search_columns
    )
    df_cleaned["set"] = df_cleaned.apply(
        _map_names, axis=1, mapping=set_mapping, search_columns=search_columns
    ).fillna(set_mapping["_"])

    df_cleaned = df_cleaned.dropna(subset=["carrier"])

    df_capacity = (
        df_cleaned.groupby(["bus", "year", "carrier", "set"])["data"]
        .sum()
        .rename("p_nom")
        .reset_index()
    )

    return df_capacity


def add_eur_capacities(
    df: pd.DataFrame,
    carrier_mapping: dict,
    set_mapping: dict,
) -> pd.DataFrame:
    """
    Function to append aggregated capacities for each European countries represented in the model.

    Args:
        df (pd.DataFrame): European country-level aggregated powerplant dataframe
        carrier_mapping (dict): dictionary to map technologies to PyPSA-Eur carriers names
        set_mapping (dict): dictionary to map technologies to powerplantmatching set names
    """
    search_columns = ["Type", "SubType", "SubSubType"]
    df_cleaned = df[df["Variable"] == "Capacity (MW)"]
    df_cleaned["carrier"] = df_cleaned.apply(
        _map_names, axis=1, mapping=carrier_mapping, search_columns=search_columns
    )

    df_cleaned["set"] = df_cleaned.apply(
        _map_names, axis=1, mapping=set_mapping, search_columns=search_columns
    ).fillna(set_mapping["_"])

    df_cleaned = df_cleaned.dropna(subset=["carrier"])
    df_capacity = (
        df_cleaned.set_index(["bus", "year", "carrier", "set"])["data"]
        .rename("p_nom")
        .reset_index()
    )
    return df_capacity


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("create_powerplants_table")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the file paths
    df_gsp = pd.read_csv(snakemake.input.gsp_data).query("Template == 'Generation'")
    df_eur = pd.read_csv(snakemake.input.eur_data)

    # Load all the params
    carrier_mapping_gb = snakemake.params.carrier_mapping_gb
    carrier_mapping_eur = snakemake.params.carrier_mapping_eur
    set_mapping = snakemake.params.set_mapping

    df_capacity_gb = table_gb_capacities(df_gsp, carrier_mapping_gb, set_mapping)
    logger.info("Tabulated the capacities into a table in PyPSA-Eur format")

    df_capacity_eur = add_eur_capacities(
        df_eur,
        carrier_mapping_eur,
        set_mapping,
    )
    logger.info("Added the EU wide capacities to the capacity table")

    df_capacity = pd.concat([df_capacity_gb, df_capacity_eur], ignore_index=True)
    df_capacity.to_csv(snakemake.output.csv, index=False)
