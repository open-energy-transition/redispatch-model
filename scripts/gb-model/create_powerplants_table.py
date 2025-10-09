# SPDX-FileCopyrightText:  gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


"""
Capacity table generator.

This is a script to extract powerplant capacities at the GSP level from the BB1 sheet of the FES workbook
"""

import pandas as pd
import numpy as np
import logging
from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def parse_inputs(bb1_path, bb2_path, fes_scenario, fes_year):
    """
    Parse the input data to the required format

    Args:
    1. bb1_path - path of extracted sheet BB1 of the FES workbook
    2. bb2_path - path of extracted sheet BB2 of the FES workbook
    3. fes_scenario - FES scenario  
    4. fes_year - Year for which model is built
    """     

    df_bb2=pd.read_csv(bb2_path)
    df_bb2=df_bb2.query("Template == 'Generation '")
    df_bb2.loc[:,'Techname']=df_bb2.apply(lambda x: x['Technology'] + " " + x['Technology Detail'], axis=1)
    df_bb2=df_bb2.loc[df_bb2['Parameter'] == 'Building Block ID Number']
    bbid_tech_map=dict(df_bb2[['data','Techname']].values)

    df_bb1=pd.read_csv(bb1_path)
    df_bb1_ltw=df_bb1.query("`FES Scenario` == @fes_scenario")\
                    .query("Unit == 'MW'")\
                    .query('`Building Block ID Number` == @keys',local_dict={'keys':list(bbid_tech_map.keys())})\
                    .query('year == @fes_year')

    return bbid_tech_map, df_bb1_ltw


def table_gb_capacities(df, bbid_tech_map):
    """
    To table the powerplant capacities in a format required by PyPSA-Eur

    Args:
    1. df - FES powerplant data
    2. bbid_tech_map - dictionary to map Building Block ID numbers to the technology in FES workbook
    """ 

    df_capacity=pd.DataFrame(index=['Name','Fueltype','Technology','Set',
                                    'Country','Capacity','Efficiency','Duration',
                                    'Volume_Mm3','DamHeight_m','StorageCapacity_MWh',
                                    'DateIn','DateRetrofit','DateMothball','DateOut',
                                    'lat','lon','EIC','projectID'])

    df = df.drop(columns=['FES Scenario','year','Unit'])
    for bb_id in bbid_tech_map.keys():
        df_ppl=df.loc[df['Building Block ID Number'] == bb_id]
        tech_detail=bbid_tech_map[bb_id]
        df_ppl['Name']=bb_id+" "+df_ppl['GSP']
        df_ppl.rename(columns={'data':'Capacity'},inplace=True)
        if "CHP" in tech_detail:
            df_ppl['Set'] = 'CHP'
        elif "hydro" in tech_detail:
            df_ppl['Set'] = 'Reservoir' #Check if this needs to be changed to run of river
        else:
            df_ppl['Set'] = 'PP'

        df_capacity=df_capacity._append(df_ppl, ignore_index=True)
    df_capacity['Country']='GB'

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("create_powerplants_table")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    bb1_path=snakemake.input.bb1_sheet
    bb2_path=snakemake.input.bb2_sheet
    gsp_coordinates_path=snakemake.input.gsp_coordinates
    fes_scenario=snakemake.params.scenario
    year=snakemake.params.year

    bbid_tech_map, df_bb1_ltw=parse_inputs(bb1_path, bb2_path, fes_scenario, year)
    logger.info(f"Extracted the {fes_scenario} relevant data")

    table_gb_capacities(df_bb1_ltw, bbid_tech_map)
    logger.info("Tabulated the capacities into a table in PyPSA-Eur format")