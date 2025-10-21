# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

"""
    Hydrogen data processor.

    This script processes hydrogen data from FES workbook.
"""

import logging
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).parent))
from scripts._helpers import configure_logging, set_scenario_config
from _helpers import _strip_str

logger = logging.getLogger(__name__)


def _standartize_year(series: pd.Series) -> pd.Series:
    """Standartize year format in a pandas Series."""
    if series.dtype == 'object' and '-' in str(series.iloc[0]):
        series = pd.to_datetime(series).dt.year
    return series.astype(int) if series.dtype == "object" else series


def parse_demand_inputs(
    demand_sheets: list,
    demand_sheets_mapping: dict,
    fes_scenario: str,
    year_range: list,
) -> pd.DataFrame:
    """
    Parse the hydrogen demand data to the required format.

    Args:
        demand_sheets (list): List of file paths to demand sheet CSV files
        demand_sheets_mapping (dict): Nested mapping of fes_year -> sector -> sheet_name
        fes_scenario (str): FES scenario name to filter by (e.g., "leading the way")
        year_range (list): Two-element list [start_year, end_year] defining the year range to include

    Returns:
        df_hydrogen_demand (pd.DataFrame): Processed hydrogen demand data with year as index
        and sectors as columns. Data values are in TWh for all sectors.
    """

    # Retrieve annual demand data from sheets
    demand_data = {}
    df_hydrogen_demand = pd.DataFrame()
    # Iterate through the nested structure: year -> sector -> sheet_name
    for year, sectors in demand_sheets_mapping.items():
        for sector, sheet_name in sectors.items():
            # Find the corresponding file in the demand_sheets
            # Look for files that contain both the year and sheet_name
            sheet_file = next(
                f for f in demand_sheets
                if sheet_name in f and str(year) in f
            )

            demand_data[sector] = pd.read_csv(sheet_file)

            # Remove whitespace from string columns
            demand_data[sector] = demand_data[sector].apply(_strip_str)

            # Make all column names lowercase
            demand_data[sector].columns = demand_data[sector].columns.str.lower()

            # Standartize year format
            demand_data[sector]['year'] = _standartize_year(demand_data[sector]['year'])

            # Make data column numeric
            demand_data[sector]['data'] = pd.to_numeric(
                demand_data[sector]['data'], errors='coerce'
            ).fillna(0)

            # Filter by year range
            demand_data[sector] = demand_data[sector][
                demand_data[sector]['year'].between(year_range[0], year_range[1])
            ]

            # Add scenario column for Road Transport
            if sector == 'road transport' and 'scenario' not in demand_data[sector].columns:
                demand_data[sector]['scenario'] = fes_scenario

            # Select scenario
            demand_data[sector] = demand_data[sector][
                demand_data[sector]['scenario'].str.lower() == fes_scenario
            ]

            # Select hydrogen demand for road transport
            if sector == 'road transport':
                demand_data[sector] = demand_data[sector][
                    demand_data[sector]['carrier'].str.lower() == 'hydrogen'
                ]

            # Convert demand to TWh
            if sector in ["industrial", "commercial"]:
                demand_data[sector]['data'] = demand_data[sector]['data'] / 1e3  # GWh to TWh

            # Define hydrogen demand dataframe
            if sector != 'other':
                df_hydrogen_demand[sector] = demand_data[sector][['year', 'data']].set_index('year')
            elif sector == 'other':
                df_sector = demand_data[sector]

                for sub_sector in other_sectors_list:
                    df_hydrogen_demand[sub_sector] = df_sector[
                        (df_sector["fuel type"].str.lower() == "hydrogen") &
                        (df_sector["category lookup"].str.lower() == sub_sector) &
                        (df_sector["category"].str.lower() == "demand")
                    ].groupby("year")["data"].sum()

    return df_hydrogen_demand


def get_regional_distribution(df: pd.Series) -> pd.Series:
    """
    Calculate regional distribution of data for each year.

    Args:
        df (pd.Series): Series containing data with index as 'bus' and 'year'.

    Returns:
        pd.Series: Series with the same index as input, but containing regional distribution
                   proportions instead of absolute values. Each row (year) sums to 1.0 across all
                   regions (columns).
    """
    # Calculate totals per year
    yearly_totals = df.groupby("year").sum()

    # Calculate regional distribution
    regional_distribution = df / yearly_totals

    return regional_distribution


def parse_grid_electrolysis_data(
    regional_gb_data_path: str,
) -> pd.DataFrame:
    """
    Parse and process grid-connected hydrogen electrolysis supply data.

    This function processes regional hydrogen electrolysis capacity data, calculates 
    regional distributions, and redistributes unmapped capacities based on existing
    regional patterns.

    Args:
        regional_gb_data_path (str): Path to the regional GB data CSV file containing
                                   hydrogen electrolysis capacity data by region and year.

    Returns:
        pd.DataFrame: Processed grid-connected electrolysis capacities with MultiIndex 
                     ['bus', 'year'] and 'data' values in original units in MW.
                     Includes both originally mapped capacities and redistributed 
                     unmapped capacities.

    Processing steps:
        1. Filter regional data for hydrogen electrolysis technology
        2. Group mapped data by bus (region) and year
        3. Calculate regional distribution proportions for each year
        4. Identify and aggregate unmapped capacities (NaN bus values)
        5. Redistribute unmapped capacities based on regional distribution patterns
        6. Combine mapped and redistributed capacities into final dataset
    """
    # Read regional GB data
    regional_gb_data = pd.read_csv(regional_gb_data_path)

    # Select regional hydrogen electrolysis data
    electrolysis_data = regional_gb_data[
        (regional_gb_data["Technology Detail"].str.lower() == "hydrogen electrolysis")
    ]

    # Calculate regional grid-connected electrolysis capacities
    regional_grid_electrolysis_capacities = electrolysis_data.groupby(
        ["bus", "year"]
    )["data"].sum()

    # Calculate regional distribution of grid-connected electrolysis capacities
    electrolysis_distribution = get_regional_distribution(
        regional_grid_electrolysis_capacities
    )

    # Map unmapped grid-connected electrolysis data
    unmapped_grid_electrolysis_capacities = electrolysis_data[
        electrolysis_data["bus"].isnull()
    ].groupby("year")["data"].sum()
    unmapped_grid_electrolysis_capacities = unmapped_grid_electrolysis_capacities * electrolysis_distribution

    # Combine mapped and unmapped grid-connected electrolysis data
    grid_electrolysis_capacities = regional_grid_electrolysis_capacities + unmapped_grid_electrolysis_capacities

    # Rename series to 'p_nom'
    grid_electrolysis_capacities.name = "p_nom"

    return grid_electrolysis_capacities


def parse_supply_data(
    supply_sheets: list,
    supply_sheets_mapping: dict,
    year_range: list,
) -> pd.DataFrame:
    """
    Parse the hydrogen supply data to the required format.

    Args:
        supply_sheets (list): List of file paths to supply sheet CSV files
        supply_sheets_mapping (dict): Nested mapping of fes_year -> data type -> sheet_name
        year_range (list): Two-element list [start_year, end_year] defining the year range to include

    Returns:
        pd.DataFrame: Processed hydrogen supply data with MultiIndex ['bus', 'year'] and 'data' values in original units in TWh.
    """
    # Retrieve annual supply data from sheets
    for year, sheets in supply_sheets_mapping.items():
        for data_type, sheet_name in sheets.items():
            # Find the corresponding file in the supply_sheets
            if data_type == "hydrogen_supply":
                sheet_file = next(
                    f for f in supply_sheets if sheet_name in f
                )

    # Read the supply sheet data
    supply_data = pd.read_csv(sheet_file)
    supply_data = supply_data.apply(_strip_str)

    # Select exogeneous hydrogen supply data
    supply_data = supply_data[
        supply_data["carrier"].str.lower().isin(
            exogeneous_supply_list
        )
    ]

    # Standartize year format
    supply_data['year'] = _standartize_year(supply_data['year'])

    # Make data column numeric
    supply_data['data'] = pd.to_numeric(
        supply_data['data'], errors='coerce'
    ).fillna(0)

    # Filter by year range
    supply_data = supply_data[
        supply_data['year'].between(year_range[0], year_range[1])
    ]

    # Use lowercase for carrier names
    supply_data['carrier'] = supply_data['carrier'].str.lower()

    # Pivot the data: year as index, carriers as columns, data as values
    df_hydrogen_supply = supply_data.pivot_table(
        index='year',
        columns='carrier',
        values='data',
        aggfunc='sum',  # Sum if there are duplicates
        fill_value=0    # Fill missing combinations with 0
    )

    return df_hydrogen_supply


def parse_electricity_demand(
    supply_sheets: list,
    supply_sheets_mapping: dict,
    year_range: list,
    grid_electrolysis_capacities: pd.DataFrame,
) -> pd.DataFrame:
    """
    Parse the electricity demand data of non-grid electrolysis.

    This function processes electricity demand data for industrial and commercial
    off-grid electrolysis and distributes it regionally based on grid-connected
    electrolysis capacity patterns.

    Args:
        supply_sheets (list): List of file paths to supply sheet CSV files
        supply_sheets_mapping (dict): Nested mapping of fes_year -> data_type -> sheet_name
        year_range (list): Two-element list [start_year, end_year] defining the year range to include
        grid_electrolysis_capacities (pd.DataFrame): Grid-connected electrolysis capacities
                                                    with MultiIndex ['bus', 'year'] for regional distribution

    Returns:
        pd.DataFrame: Processed electricity demand data with MultiIndex ['bus', 'year']
                     and electricity demand values in TWh. Data is regionally distributed
                     based on grid-connected electrolysis capacity patterns.

    Processing steps:
        1. Load electricity demand data from supply sheets
        2. Filter for "I&C off grid electrolysis demand" data items
        3. Standardize year format and convert data to numeric
        4. Filter by specified year range
        5. Convert units from GWh to TWh
        6. Calculate regional distribution based on grid electrolysis capacities
        7. Apply regional distribution to electricity demand data
    """
    # Retrieve annual electricity demand for non-grid electrolysis data
    for year, sheets in supply_sheets_mapping.items():
        for data_type, sheet_name in sheets.items():
            # Find the corresponding file in the supply_sheets
            if data_type == "electricity_demand":
                sheet_file = next(
                    f for f in supply_sheets if sheet_name in f
                )

    # Read the electricity demand sheet data
    electricity_demand_data = pd.read_csv(sheet_file)
    electricity_demand_data = electricity_demand_data.apply(_strip_str)

    # Make columns lowercase
    electricity_demand_data.columns = electricity_demand_data.columns.str.lower()

    # Standardize year format
    electricity_demand_data['year'] = _standartize_year(electricity_demand_data['year'])

    # Make data column numeric
    electricity_demand_data['data'] = pd.to_numeric(
        electricity_demand_data['data'], errors='coerce'
    ).fillna(0)

    # Filter by year range
    electricity_demand_data = electricity_demand_data[
        electricity_demand_data['year'].between(year_range[0], year_range[1])
    ]

    # Select non-grid electrolysis electricity demand
    electricity_demand_data = electricity_demand_data[
        electricity_demand_data["data item"].str.lower() == "i&c off grid electrolysis demand (lw only)"
    ]

    # Prepare final dataframe
    df_electricity_demand = electricity_demand_data.set_index('year')["data"]

    # Convert to TWh
    df_electricity_demand = df_electricity_demand / 1e3  # GWh to TWh

    # Calculate regional distribution of grid-connected electrolysis capacities
    electrolysis_distribution = get_regional_distribution(
        grid_electrolysis_capacities
    )

    # Apply regional distribution to electricity demand data
    df_electricity_demand = df_electricity_demand * electrolysis_distribution

    return df_electricity_demand


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the file path
    demand_sheets = snakemake.input.demand_sheets
    regional_gb_data_path = snakemake.input.regional_gb_data
    supply_sheets = snakemake.input.supply_sheets

    # Load all params
    fes_scenario = snakemake.params.scenario
    year_range = snakemake.params.year_range
    demand_sheets_mapping = snakemake.params.fes_demand_sheets
    other_sectors_list = snakemake.params.other_sectors_list
    supply_sheets_mapping = snakemake.params.fes_supply_sheets
    exogeneous_supply_list = snakemake.params.exogeneous_supply_list

    # Parse demand data
    df_hydrogen_demand = parse_demand_inputs(
        demand_sheets,
        demand_sheets_mapping,
        fes_scenario,
        year_range,
    )

    # Save the hydrogen demand data
    df_hydrogen_demand.to_csv(snakemake.output.hydrogen_demand)

    # Parse grid connected electrolysis data
    grid_electrolysis_capacities = parse_grid_electrolysis_data(
        regional_gb_data_path,
    )

    # Save the grid-connected electrolysis capacities
    grid_electrolysis_capacities.to_csv(
        snakemake.output.grid_electrolysis_capacities
    )

    # Parse hydrogen supply data
    df_hydrogen_supply = parse_supply_data(
        supply_sheets,
        supply_sheets_mapping,
        year_range,
    )

    # Save the hydrogen supply data
    df_hydrogen_supply.to_csv(snakemake.output.hydrogen_supply)

    # Parse electricity demand data for non-grid electrolysis
    df_electricity_demand = parse_electricity_demand(
        supply_sheets,
        supply_sheets_mapping,
        year_range,
        grid_electrolysis_capacities,
    )

    # Save the electricity demand data
    df_electricity_demand.to_csv(snakemake.output.electricity_demand)
