# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


import logging

import geopandas as gpd
import numpy as np
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
        geometry=gpd.points_from_xy(df[lon_col], df[lat_col]), crs=crs
    ).to_crs(gdf_regions.crs)
    regions = gpd.sjoin(points, gdf_regions, how="left", predicate="intersects")["name"]
    return regions


def strip_srt(series: pd.Series) -> pd.Series:
    """Strip whitespace from strings in a pandas Series."""
    return series.str.strip() if series.dtype == "object" else series


def to_numeric(series: pd.Series) -> pd.Series:
    """Convert a pandas Series to numeric, replacing - with 0."""
    series = series.astype(str).str.strip()
    series = series.replace("-", np.nan)
    series = pd.to_numeric(series).fillna(0)
    return series
