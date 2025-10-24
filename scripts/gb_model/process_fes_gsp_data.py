# SPDX-FileCopyrightText: gb-dispatch-model contributors
#
# SPDX-License-Identifier: MIT


"""
GSP-level data table generator.

This is a script to combine the BB1 sheet with the BB2 (metadata) sheet of the FES workbook.
"""

import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import map_points_to_regions

logger = logging.getLogger(__name__)


def _strip_str(series: pd.Series) -> pd.Series:
    """Strip whitespace from strings in a pandas Series."""
    return series.str.strip() if series.dtype == "object" else series


def parse_inputs(
    bb1_path: str,
    bb2_path: str,
    gsp_coordinates_path: str,
    fes_scenario: str,
    year_range: list,
) -> pd.DataFrame:
    """
    Parse the input data to the required format.

    Args:
        bb1_path (str): path of extracted sheet BB1 of the FES workbook
        bb2_path (str): path of extracted sheet BB2 of the FES workbook
        gsp_coordinates_path (str): path of the GSP supply point coordinates file
        fes_scenario (str): FES scenario
    """

    df_bb2 = pd.read_csv(bb2_path)

    # First step: extract the ID numbers from the Parameter column and set it as the index (it is the only unique identifier for table BB2)
    df_bb2 = (
        df_bb2.set_index(
            ["Template", "Technology", "Technology Detail", "Parameter"], append=True
        )
        .squeeze()
        .unstack("Parameter")
    )
    df_bb2_pivoted = (
        df_bb2.bfill()
        .where(~df_bb2["Building Block ID Number"].isnull())
        .dropna(how="all")
        .reset_index()
        .set_index("Building Block ID Number")
        .drop("level_0", axis=1)
        .apply(_strip_str)
    )

    df_bb1 = pd.read_csv(bb1_path)
    df_bb1 = df_bb1.apply(_strip_str)
    df_bb1_scenario = df_bb1[
        (df_bb1["FES Scenario"] == fes_scenario)
        & (df_bb1["year"].isin(range(year_range[0], year_range[1] + 1)))
    ]
    df_bb1_bb2_scenario = pd.merge(
        df_bb1_scenario,
        df_bb2_pivoted,
        left_on="Building Block ID Number",
        right_index=True,
    )
    assert len(df_bb1_bb2_scenario) == len(df_bb1_scenario), (
        "Some Building Blocks in BB1 are not present in BB2"
    )

    # We allow cases where there is only a partial match ("Number" vs "Number of" by comparing string starts)
    units_match = df_bb1_bb2_scenario.apply(
        lambda x: x.Units.startswith(x.Unit), axis=1
    )
    assert (units_match).all(), (
        "Mapping of building blocks between BB1 and BB2 may be incorrect as some units do not match: "
        f"{df_bb1_bb2_scenario[~units_match][['Unit', 'Units']]}"
    )

    df_bb1_bb2_scenario = df_bb1_bb2_scenario.drop(
        columns=["Units", "Building Block ID Number"]
    )

    df_gsp_coordinates = pd.read_csv(gsp_coordinates_path)
    df_gsp_coordinates = df_gsp_coordinates.apply(_strip_str)

    # Note
    # The GSP's "East Claydon" and "Ferrybridge B" have duplicates
    # the lat and lon information is the same but the GSP ID and GSP group are slightly different
    df_gsp_coordinates = df_gsp_coordinates.drop_duplicates(subset=["Name"])
    df_bb1_bb2_with_lat_lon = pd.merge(
        df_bb1_bb2_scenario, df_gsp_coordinates, left_on="GSP", right_on="Name"
    )

    # Missing data checks.
    # We won't raise errors here as we are willing to accept some missing data for now
    missing_lat_lon = df_bb1_bb2_with_lat_lon[
        df_bb1_bb2_with_lat_lon[["Latitude", "Longitude"]].isnull().any(axis=1)
    ].GSP.unique()
    if len(missing_lat_lon) > 0:
        logger.error(
            f"The following GSPs are missing latitude and/or longitude information: {missing_lat_lon}"
        )

    missing_gsps = set(df_bb1_bb2_scenario.GSP).difference(df_bb1_bb2_with_lat_lon.GSP)
    if missing_gsps:
        logger.error(
            f"The following GSPs are missing from the GSP coordinates file: {missing_gsps}"
        )
    df_final = pd.concat(
        [df_bb1_bb2_scenario.query("GSP in @missing_gsps"), df_bb1_bb2_with_lat_lon]
    )
    return df_final


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the file paths
    bb1_path = snakemake.input.bb1_sheet
    bb2_path = snakemake.input.bb2_sheet
    gsp_coordinates_path = snakemake.input.gsp_coordinates
    gdf_regions = gpd.read_file(snakemake.input.regions)

    # Load all the params
    fes_scenario = snakemake.params.scenario
    year_range = snakemake.params.year_range

    df = parse_inputs(
        bb1_path, bb2_path, gsp_coordinates_path, fes_scenario, year_range
    )
    df["bus"] = map_points_to_regions(df, gdf_regions, "Latitude", "Longitude")
    logger.info(f"Extracted the {fes_scenario} relevant data")

    df.to_csv(snakemake.output.csv, index=False)
