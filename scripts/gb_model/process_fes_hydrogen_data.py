# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

"""
Hydrogen data processor.

This script processes hydrogen data from FES workbook.
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from scripts._helpers import configure_logging, set_scenario_config
from scripts.gb_model._helpers import strip_srt, to_numeric

logger = logging.getLogger(__name__)


def _standardize_year(series: pd.Series) -> pd.Series:
    """Standardize year format in a pandas Series."""
    if series.dtype == "object" and "-" in str(series.iloc[0]):
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
        and sectors as columns. Data values are in MWh for all sectors.
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
                f for f in demand_sheets if sheet_name in f and str(year) in f
            )

            demand_data[sector] = pd.read_csv(sheet_file)

            # Remove whitespace from string columns
            demand_data[sector] = demand_data[sector].apply(strip_srt)

            # Make all column names lowercase
            demand_data[sector].columns = demand_data[sector].columns.str.lower()

            # Standardize year format
            demand_data[sector]["year"] = _standardize_year(demand_data[sector]["year"])

            # Make data column numeric
            demand_data[sector]["data"] = to_numeric(demand_data[sector]["data"])

            # Filter by year range
            demand_data[sector] = demand_data[sector][
                demand_data[sector]["year"].between(year_range[0], year_range[1])
            ]

            # Add scenario column for Road Transport
            if "scenario" not in demand_data[sector].columns:
                demand_data[sector]["scenario"] = fes_scenario

            # Select scenario
            demand_data[sector] = demand_data[sector][
                demand_data[sector]["scenario"].str.lower() == fes_scenario
            ]

            # Select hydrogen demand for road transport
            if sector == "road transport":
                demand_data[sector] = demand_data[sector][
                    demand_data[sector]["carrier"].str.lower() == "hydrogen"
                ]

            # Convert demand to MWh
            if sector in ["industrial", "commercial"]:
                demand_data[sector]["data"] = (
                    demand_data[sector]["data"] * 1e3
                )  # GWh to MWh
            else:
                demand_data[sector]["data"] = (
                    demand_data[sector]["data"] * 1e6
                )  # TWh to MWh

            # Define hydrogen demand dataframe
            if sector != "other":
                df_hydrogen_demand[sector] = demand_data[sector][
                    ["year", "data"]
                ].set_index("year")
            elif sector == "other":
                df_sector = demand_data[sector]

                for sub_sector in other_sectors_list:
                    df_hydrogen_demand[sub_sector] = (
                        df_sector[
                            (df_sector["fuel type"].str.lower() == "hydrogen")
                            & (df_sector["category lookup"].str.lower() == sub_sector)
                            & (df_sector["category"].str.lower() == "demand")
                        ]
                        .groupby("year")["data"]
                        .sum()
                    )

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
    regional_distribution = df.groupby("year").apply(lambda x: x / x.sum())

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
    regional_grid_electrolysis_capacities = electrolysis_data.groupby(["bus", "year"])[
        "data"
    ].sum()

    # Calculate regional distribution of grid-connected electrolysis capacities
    electrolysis_distribution = get_regional_distribution(
        regional_grid_electrolysis_capacities
    )

    # Map unmapped grid-connected electrolysis data
    unmapped_grid_electrolysis_capacities = (
        electrolysis_data[electrolysis_data["bus"].isnull()]
        .groupby("year")["data"]
        .sum()
    )
    unmapped_grid_electrolysis_capacities = (
        unmapped_grid_electrolysis_capacities * electrolysis_distribution
    )

    # Combine mapped and unmapped grid-connected electrolysis data
    grid_electrolysis_capacities = (
        regional_grid_electrolysis_capacities + unmapped_grid_electrolysis_capacities
    )

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
        pd.DataFrame: Processed hydrogen supply data with MultiIndex ['bus', 'year'] and 'data' values in MWh.
    """
    # Retrieve annual supply data from sheets
    for year, sheets in supply_sheets_mapping.items():
        for data_type, sheet_name in sheets.items():
            # Find the corresponding file in the supply_sheets
            if data_type == "hydrogen_supply":
                sheet_file = next(f for f in supply_sheets if sheet_name in f)

    # Read the supply sheet data
    supply_data = pd.read_csv(sheet_file)
    supply_data = supply_data.apply(strip_srt)

    # Select exogeneous hydrogen supply data
    supply_data = supply_data[
        supply_data["carrier"].str.lower().isin(exogeneous_supply_list)
    ]

    # Standardize year format
    supply_data["year"] = _standardize_year(supply_data["year"])

    # Make data column numeric
    supply_data["data"] = to_numeric(supply_data["data"])

    # Filter by year range
    supply_data = supply_data[supply_data["year"].between(year_range[0], year_range[1])]

    # Use lowercase for carrier names
    supply_data["carrier"] = supply_data["carrier"].str.lower()

    # Pivot the data: year as index, carriers as columns, data as values
    df_hydrogen_supply = supply_data.pivot_table(
        index="year",
        columns="carrier",
        values="data",
        aggfunc="sum",  # Sum if there are duplicates
        fill_value=0,  # Fill missing combinations with 0
    )

    # Convert to MWh
    df_hydrogen_supply = df_hydrogen_supply * 1e6  # TWh to MWh

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
                     and electricity demand values in MWh. Data is regionally distributed
                     based on grid-connected electrolysis capacity patterns.

    Processing steps:
        1. Load electricity demand data from supply sheets
        2. Filter for "I&C off grid electrolysis demand" data items
        3. Standardize year format and convert data to numeric
        4. Filter by specified year range
        5. Convert units from GWh to MWh
        6. Calculate regional distribution based on grid electrolysis capacities
        7. Apply regional distribution to electricity demand data
    """
    # Retrieve annual electricity demand for non-grid electrolysis data
    for year, sheets in supply_sheets_mapping.items():
        for data_type, sheet_name in sheets.items():
            # Find the corresponding file in the supply_sheets
            if data_type == "electricity_demand":
                sheet_file = next(f for f in supply_sheets if sheet_name in f)

    # Read the electricity demand sheet data
    electricity_demand_data = pd.read_csv(sheet_file)
    electricity_demand_data = electricity_demand_data.apply(strip_srt)

    # Make columns lowercase
    electricity_demand_data.columns = electricity_demand_data.columns.str.lower()

    # Standardize year format
    electricity_demand_data["year"] = _standardize_year(electricity_demand_data["year"])

    # Make data column numeric
    electricity_demand_data["data"] = to_numeric(electricity_demand_data["data"])

    # Filter by year range
    electricity_demand_data = electricity_demand_data[
        electricity_demand_data["year"].between(year_range[0], year_range[1])
    ]

    # Select non-grid electrolysis electricity demand
    electricity_demand_data = electricity_demand_data[
        electricity_demand_data["data item"].str.lower()
        == "i&c off grid electrolysis demand (lw only)"
    ]

    # Prepare final dataframe
    df_electricity_demand = electricity_demand_data.set_index("year")["data"]

    # Convert to MWh
    df_electricity_demand = df_electricity_demand * 1e3  # GWh to MWh

    # Calculate regional distribution of grid-connected electrolysis capacities
    electrolysis_distribution = get_regional_distribution(grid_electrolysis_capacities)

    # Apply regional distribution to electricity demand data
    df_electricity_demand = df_electricity_demand * electrolysis_distribution

    # Rename series to 'p_set'
    df_electricity_demand.name = "p_set"

    return df_electricity_demand


def parse_storage_data(
    storage_sheets: list,
    storage_sheets_mapping: dict,
    fes_scenario: str,
    year_range: list,
    interpolation_method: str = "linear",
) -> pd.DataFrame:
    """
    Parse the hydrogen storage data to the required format.

    Args:
        storage_sheets (list): List of file paths to storage sheet CSV files
        storage_sheets_mapping (dict): Nested mapping of fes_year -> data type -> sheet_name
        fes_scenario (str): FES scenario name to filter by (e.g., "Leading the Way")
        year_range (list): Two-element list [start_year, end_year] defining the year range to include
        interpolation_method (str): Method for interpolating between 5-year data points
                                   ("linear", "s_curve", or "step"). Default is "linear".

    Returns:
        pd.DataFrame: Processed hydrogen storage data with year as index and storage types as columns.
                     Data values are in MWh. Data is interpolated from
                     5-year intervals to annual values using the specified interpolation method.

    Processing steps:
        1. Load storage data from CSV files based on mapping configuration
        2. Filter data by specified FES scenario
        3. Standardize year format and convert data to numeric
        4. Apply interpolation to fill annual data between 5-year intervals
        5. Filter by specified year range
        6. Return processed DataFrame with annual storage capacity data e_nom in MWh
    """
    # Retrieve annual storage data from sheets
    for year, sheets in storage_sheets_mapping.items():
        for data_type, sheet_name in sheets.items():
            # Find the corresponding file in the storage_sheets
            if data_type == "storage_capacity":
                sheet_file = next(f for f in storage_sheets if sheet_name in f)

    # Read the storage sheet data
    storage_data = pd.read_csv(sheet_file)
    storage_data = storage_data.apply(strip_srt)

    # Select the scenario data
    storage_data = storage_data[storage_data["scenario"].str.lower() == fes_scenario]

    # Standardize year format
    storage_data["year"] = _standardize_year(storage_data["year"])

    # Make data column numeric
    storage_data["data"] = to_numeric(storage_data["data"])

    # Perform interpolation to get annual data
    storage_data_interpolated = interpolate_yearly_data(
        storage_data,
        method=interpolation_method,
        year_column="year",
        data_column="data",
    )

    # Filter by year range
    storage_data_interpolated = storage_data_interpolated[
        storage_data_interpolated["year"].between(year_range[0], year_range[1])
    ]

    # Prepare final dataframe
    storage_data_interpolated = storage_data_interpolated.set_index("year")["data"]
    storage_data_interpolated.name = "e_nom"

    return storage_data_interpolated


def interpolate_yearly_data(
    df: pd.DataFrame,
    method: str = "linear",
    year_column: str = "year",
    data_column: str = "data",
) -> pd.DataFrame:
    """
    Interpolate data between 5-year intervals to get annual values.

    Args:
        df (pd.DataFrame): DataFrame with sparse year data (e.g., every 5 years)
        method (str): Interpolation method - "linear", "s_curve", or "step"
        year_column (str): Name of the year column
        data_column (str): Name of the data column to interpolate

    Returns:
        pd.DataFrame: DataFrame with interpolated annual data
    """

    # Create full year range
    min_year = df[year_column].min()
    max_year = df[year_column].max()
    full_years = np.arange(min_year, max_year + 1)

    # Get the sparse data points
    years = df[year_column].values
    values = df[data_column].values

    if method == "linear":
        # Linear interpolation
        f = interp1d(years, values, kind="linear", fill_value="extrapolate")
        interpolated_values = f(full_years)

    elif method == "s_curve":
        # S-curve (sigmoid-like) interpolation using cubic spline
        f = interp1d(years, values, kind="cubic", fill_value="extrapolate")
        interpolated_values = f(full_years)

        # Ensure non-negative values (optional)
        interpolated_values = np.maximum(interpolated_values, 0)

    elif method == "step":
        # Step interpolation - hold previous value until next data point
        f = interp1d(years, values, kind="previous", fill_value="extrapolate")
        interpolated_values = f(full_years)

    else:
        raise ValueError("Method must be 'linear', 's_curve', or 'step'")

    # Create new DataFrame with interpolated data
    result_df = pd.DataFrame(
        {year_column: full_years, data_column: interpolated_values}
    )

    # Add other columns if they exist (like scenario)
    for col in df.columns:
        if col not in [year_column, data_column]:
            # Use the first value for other columns (assuming they're constant)
            result_df[col] = df[col].iloc[0]

    return result_df


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
    storage_sheets = snakemake.input.storage_sheet

    # Load all params
    fes_scenario = snakemake.params.scenario
    year_range = snakemake.params.year_range
    demand_sheets_mapping = snakemake.params.fes_demand_sheets
    other_sectors_list = snakemake.params.other_sectors_list
    supply_sheets_mapping = snakemake.params.fes_supply_sheets
    exogeneous_supply_list = snakemake.params.exogeneous_supply_list
    storage_sheets_mapping = snakemake.params.fes_storage_sheets
    interpolation_method = snakemake.params.interpolation_method

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
    grid_electrolysis_capacities.to_csv(snakemake.output.grid_electrolysis_capacities)

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

    # Parse storage data
    df_hydrogen_storage = parse_storage_data(
        storage_sheets,
        storage_sheets_mapping,
        fes_scenario,
        year_range,
        interpolation_method,
    )

    # Save the storage data
    df_hydrogen_storage.to_csv(snakemake.output.hydrogen_storage)
