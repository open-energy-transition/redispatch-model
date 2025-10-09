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
import re

logger = logging.getLogger(__name__)

# Rename fueltype to PyPSA-Eur convention
fuel_dict={
    'Biomass & Energy Crops (including CHP) ':'biomass',
    'CCGTs (non CHP)':'CCGT',
    'Coal':'coal',
    'Geothermal':'geothermal',
    'Hydro':'hydro',
    'Nuclear':'nuclear',
    'OCGTs (non CHP)':'OCGT',
    'Solar Generation':'solar'
}

def parse_inputs(bb1_path, bb2_path, gsp_coordinates_path, fes_scenario, fes_year):
    """
    Parse the input data to the required format

    Args:
    * bb1_path - path of extracted sheet BB1 of the FES workbook
    * bb2_path - path of extracted sheet BB2 of the FES workbook
    * gsp_coordinates_path - path of the GSP supply point coordinates file
    * fes_scenario - FES scenario  
    * fes_year - Year for which model is built
    """     

    df_bb2=pd.read_csv(bb2_path)
    df_bb2=df_bb2.query("Template == 'Generation '")
    df_bb2=df_bb2.loc[df_bb2['Parameter'] == 'Building Block ID Number']
    df_bb2[['Technology','Technology Detail']]=df_bb2[['Technology','Technology Detail']].apply(lambda x: x.str.strip(), axis=1)
    bbid_tech_map=df_bb2.set_index('data')[['Technology','Technology Detail']].to_dict(orient='index')

    df_bb1=pd.read_csv(bb1_path)
    df_bb1_ltw=df_bb1.query("`FES Scenario` == @fes_scenario")\
                    .query("Unit == 'MW'")\
                    .query('`Building Block ID Number` == @keys',local_dict={'keys':list(bbid_tech_map.keys())})\
                    .query('year == @fes_year')

    df_gsp_coordinates=pd.read_csv(gsp_coordinates_path)
    # Note
    # The GSP's "East Claydon" and "Ferrybridge B" have duplicates
    # the lat and lon information is the same but the GSP ID and GSP group are slightly different
    df_gsp_coordinates=df_gsp_coordinates.drop_duplicates(subset=['Name'])
    gsp_coord_dict=df_gsp_coordinates.set_index('Name')[['Latitude','Longitude']].to_dict(orient='index')

    return bbid_tech_map, df_bb1_ltw, gsp_coord_dict


def table_gb_capacities(df, bbid_tech_map, gsp_coord_dict):
    """
    To table the powerplant capacities in a format required by PyPSA-Eur

    Args:
    1. df - FES powerplant data
    2. bbid_tech_map - dictionary to map Building Block ID numbers to the technology in FES workbook
    3. gsp_coord_dict - nested dictionary of lat and lon coordinates for each GSP
    """ 

    df_capacity=pd.DataFrame(columns=['Name','Fueltype','Technology','Set',
                                    'Country','Capacity','Efficiency','Duration',
                                    'Volume_Mm3','DamHeight_m','StorageCapacity_MWh',
                                    'DateIn','DateRetrofit','DateMothball','DateOut',
                                    'lat','lon','EIC','projectID'])

    df = df.drop(columns=['FES Scenario','year','Unit'])
    for bb_id in bbid_tech_map.keys():
        df_ppl=df.loc[df['Building Block ID Number'] == bb_id]

        technology=bbid_tech_map.get(bb_id,{}).get('Technology')
        technology_detail=bbid_tech_map.get(bb_id,{}).get('Technology Detail')

        df_ppl['Name']=bb_id+" "+df_ppl['GSP']

        df_ppl['Fueltype']=fuel_dict[technology] if technology in fuel_dict.keys() else technology
        df_ppl['Technology']=technology_detail
        df_ppl.rename(columns={'data':'Capacity'},inplace=True)

        if 'non CHP' not in technology and 'CHP' in technology:
            df_ppl['Set'] = 'CHP'
        elif "hydro" in technology:
            df_ppl['Set'] = 'Reservoir' #Check if this needs to be changed to run of river
        else:
            df_ppl['Set'] = 'PP'

        df_ppl['lat']=df_ppl['GSP'].map(lambda x: gsp_coord_dict.get(x, {}).get('Latitude'))
        df_ppl['lon']=df_ppl['GSP'].map(lambda x: gsp_coord_dict.get(x, {}).get('Longitude'))

        df_ppl=df_ppl.query("Capacity != 0").query("GSP != 'Not Connected'")
        if not df_ppl.empty:
            intersecting_columns=df_ppl.columns.intersection(df_capacity.columns)
            df_capacity=pd.concat([df_capacity,df_ppl[intersecting_columns]],ignore_index=True)
    df_capacity['Country']='GB'
    df_capacity.loc[(df_capacity['Fueltype']=='Wind')&(df_capacity['Technology']=='Offshore Wind'),'Fueltype']='offwind-ac'
    df_capacity.loc[(df_capacity['Fueltype']=='Wind')&(df_capacity['Technology']=='Onshore Wind <1MW'),'Fueltype']='onwind'
    df_capacity.loc[(df_capacity['Fueltype']=='Wind')&(df_capacity['Technology']=='Onshore Wind >=1MW'),'Fueltype']='onwind'

    return df_capacity

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

    bbid_tech_map, df_bb1_ltw,gsp_coord_dict=parse_inputs(bb1_path, bb2_path, gsp_coordinates_path, fes_scenario, year)
    logger.info(f"Extracted the {fes_scenario} relevant data")

    df_capacity=table_gb_capacities(df_bb1_ltw, bbid_tech_map, gsp_coord_dict)
    logger.info("Tabulated the capacities into a table in PyPSA-Eur format")

    df_capacity.to_csv(snakemake.output.csv)