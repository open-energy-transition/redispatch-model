# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

import logging

import geopandas as gpd
import pandas as pd

logger = logging.getLogger(__name__)


def map_points_to_regions(
    df: pd.DataFrame,
    gdf_regions: gpd.GeoDataFrame,
    lat_col: str = "lat",
    lon_col: str = "lon",
    crs: str = "EPSG:4326",
) -> pd.DataFrame:
    points = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(df[lon_col], df[lat_col], index=df.index), crs=crs
    ).to_crs(gdf_regions.crs)
    regions = gpd.sjoin(points, gdf_regions, how="left", predicate="intersects").drop(
        columns="geometry"
    )
    return regions
