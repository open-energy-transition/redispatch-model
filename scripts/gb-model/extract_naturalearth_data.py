# SPDX-FileCopyrightText:  gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


"""
Extract naturalearth data.

This is a script to extract world shape files to get lat, lon information of EU countries
"""

import requests
import zipfile
import io
import geopandas as gpd
import os
import shutil
import logging
from scripts._helpers import configure_logging, set_scenario_config


logger = logging.getLogger(__name__)

def extract_files(zip_file, shape_files):
    """
        Function to extract the specific shape files from the zip archive

        Args:
            * zip_file - zip file path
            * shape_files - list of shape files to be extracted
    """
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        for file_path in shape_files:
            file=f"110m_cultural/{os.path.basename(file_path)}"
            try:
                if file in zip_ref.namelist():
                    zip_ref.extract(file, os.path.dirname(file_path))
                    # Move data one folder up
                    shutil.move(f"{os.path.dirname(shape_files[0])}/{file}",file_path)
                    logger.info(f"Successfully extracted: {file}")

                else:
                    logger.error(f"File not found in archive: {file}")
                    raise ValueError
            except Exception as e:
                logger.info(f"Error extracting {file}: {e}")

    os.rmdir(f"{os.path.dirname(shape_files[0])}/110m_cultural")

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("extract_naturalearth_data")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    zip_file=snakemake.input.zip
    shape_files=snakemake.output.shape_files

    extract_files(zip_file, shape_files)

    gdf = gpd.read_file(shape_files[0])
    gdf = gdf.to_crs(epsg=4326)

    gdf['rep_point'] = gdf.representative_point()
    output_df = gdf[['ADMIN', 'ISO_A2', 'rep_point']].copy()

    output_df['latitude'] = output_df['rep_point'].y
    output_df['longitude'] = output_df['rep_point'].x

    output_df.to_csv(snakemake.output.country_coordinates)
    logger.info("Country wise representative points extracted")



    