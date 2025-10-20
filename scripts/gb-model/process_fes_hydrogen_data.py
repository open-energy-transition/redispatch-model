# SPDX-FileCopyrightText:  gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

"""
    Hydrogen data processor.

    This script processes hydrogen data from FES workbook.
"""

import logging
from pathlib import Path

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def _standartize_year(series: pd.Series) -> pd.Series:
    """Standartize year format in a pandas Series."""
    if series.dtype == 'object' and '-' in str(series.iloc[0]):
        series = pd.to_datetime(series).dt.year
    return series.astype(int) if series.dtype == "object" else series


def parse_demand_inputs(
    bb1_path: str,
    demand_sheets: list,
    demand_sheets_mapping: dict,
    fes_scenario: str,
    year_range: list,
) -> pd.DataFrame:
    """
    Parse the hydrogen demand data to the required format.

    Args:
        hydrogen_data_path (str): path of extracted hydrogen data sheet of the FES workbook
        fes_scenario (str): FES scenario
    """

    # Retrieve annual demand data from sheets
    demand_data = {}
    for sector, sheet_name in demand_sheets_mapping.items():
        # Find the corresponding file in the demand_sheets
        sheet_file = next(f for f in demand_sheets if sheet_name in f)
        demand_data[sector] = pd.read_csv(sheet_file)

        # standartize year format
        demand_data[sector]['year'] = _standartize_year(demand_data[sector]['year'])


    # Process the hydrogen data as needed
    df_hydrogen_processed = df_hydrogen[df_hydrogen['Scenario'] == fes_scenario]

    return df_hydrogen_processed


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(Path(__file__).stem)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load the file path
    bb1_path = snakemake.input.bb1_sheet
    demand_sheets = snakemake.input.demand_sheets

    # Load all params
    fes_scenario = snakemake.params.scenario
    year_range = snakemake.params.year_range
    demand_sheets_mapping = snakemake.params.fes_demand_sheets

    # Parse demand data
    df_demand = parse_demand_inputs(
        bb1_path,
        demand_sheets,
        demand_sheets_mapping,
        fes_scenario,
        year_range,
    )
