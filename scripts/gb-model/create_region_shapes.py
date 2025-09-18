# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later


"""
Script to create region shapes by dividing country shapes based on ETYS boundary lines.

This script reads country shapes from a GeoJSON file and divides them using
boundary lines from the ETYS boundary data to create regional divisions.
The resulting regions are saved as a new GeoJSON file.
"""

import logging
import os
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.ops import split
from shapely.geometry import Point

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def load_country_shapes(filepath):
    """
    Load country shapes from GeoJSON file.

    Args:
        filepath (str): Path to the country shapes GeoJSON file

    Returns:
        geopandas.GeoDataFrame: Country shapes data
    """
    try:
        logger.info(f"Loading country shapes from: {filepath}")
        country_gdf = gpd.read_file(filepath)
        logger.info(f"Loaded {len(country_gdf)} country shapes")
        return country_gdf
    except Exception as e:
        logger.error(f"Error loading country shapes: {e}")
        raise


def load_boundary_lines(filepath):
    """
    Load boundary lines from shapefile or GeoJSON.

    Args:
        filepath (str): Path to the boundary lines file

    Returns:
        geopandas.GeoDataFrame: Boundary lines data
    """
    try:
        logger.info(f"Loading boundary lines from: {filepath}")
        boundary_gdf = gpd.read_file(filepath)
        logger.info(f"Loaded {len(boundary_gdf)} boundary features")
        return boundary_gdf
    except Exception as e:
        logger.error(f"Error loading boundary lines: {e}")
        raise


def ensure_same_crs(gdf1, gdf2):
    """
    Ensure both GeoDataFrames have the same CRS.

    Args:
        gdf1, gdf2 (geopandas.GeoDataFrame): GeoDataFrames to align

    Returns:
        tuple: Both GeoDataFrames with the same CRS
    """
    if gdf1.crs != gdf2.crs:
        logger.info(f"Converting CRS from {gdf2.crs} to {gdf1.crs}")
        gdf2 = gdf2.to_crs(gdf1.crs)
    return gdf1, gdf2


def create_regions_from_boundaries(country_shapes, boundary_lines):
    """
    Create regions by dividing country shapes using boundary lines.

    Args:
        country_shapes (geopandas.GeoDataFrame): Country polygons
        boundary_lines (geopandas.GeoDataFrame): Boundary lines for division

    Returns:
        geopandas.GeoDataFrame: Regional divisions
    """
    logger.info("Creating regions from boundaries...")

    # Convert to a projected CRS for accurate area calculations
    # Use British National Grid (EPSG:27700) which is appropriate for UK
    target_crs = "EPSG:27700"
    logger.info(f"Converting to projected CRS {target_crs} for accurate measurements")

    if country_shapes.crs != target_crs:
        country_shapes = country_shapes.to_crs(target_crs)
        logger.info(
            f"Country shapes converted. New total area: {country_shapes.geometry.area.sum() / 1000000:.0f} km²"
        )
    if boundary_lines.crs != target_crs:
        boundary_lines = boundary_lines.to_crs(target_crs)
        logger.info(
            f"Boundary lines converted. Length range: {boundary_lines.geometry.length.min():.0f} - {boundary_lines.geometry.length.max():.0f} meters"
        )

    if len(boundary_lines) == 0:
        logger.warning("No valid boundary lines found!")
        return country_shapes.copy()

    logger.info(f"Using {len(boundary_lines)} boundary lines for splitting")

    regions = []
    region_id = 1

    for idx, country in country_shapes.iterrows():
        logger.info(f"Processing country/region {idx + 1}/{len(country_shapes)}")

        try:
            # Start with the original country geometry
            current_polygons = [country.geometry]
            split_count = 0

            # Process each boundary line
            for boundary_idx, boundary in boundary_lines.iterrows():
                boundary_geom = boundary.geometry
                new_polygons = []

                for polygon in current_polygons:
                    # Check if boundary intersects this polygon
                    if polygon.intersects(boundary_geom):
                        try:
                            # Attempt to split the polygon
                            split_result = split(polygon, boundary_geom)

                            # Check if splitting actually occurred
                            if hasattr(split_result, "geoms"):
                                split_geoms = list(split_result.geoms)
                                logger.info(
                                    f"Split produced {len(split_geoms)} geometries"
                                )
                                if len(split_geoms) > 1:
                                    # Successfully split - add all valid pieces
                                    valid_pieces = 0
                                    for geom in split_geoms:
                                        if (
                                            geom.is_valid
                                            and not geom.is_empty
                                            and geom.geom_type
                                            in ["Polygon", "MultiPolygon"]
                                            and geom.area > 1000000
                                        ):  # 1 km² minimum (1,000,000 sq meters in projected CRS)
                                            new_polygons.append(geom)
                                            valid_pieces += 1
                                        else:
                                            logger.debug(
                                                f"Filtered out small geometry: area={geom.area:.0f} sq meters"
                                            )

                                    if valid_pieces > 0:
                                        split_count += 1
                                        logger.info(
                                            f"Split polygon into {valid_pieces} valid pieces (from {len(split_geoms)} total)"
                                        )
                                    else:
                                        logger.warning(
                                            f"All {len(split_geoms)} split pieces were too small - keeping original"
                                        )
                                        new_polygons.append(polygon)
                                else:
                                    # No actual split occurred
                                    new_polygons.append(polygon)
                            else:
                                # Single geometry result - no split
                                new_polygons.append(polygon)

                        except Exception as split_error:
                            logger.warning(
                                f"Split failed for boundary {boundary_idx}: {split_error}"
                            )
                            new_polygons.append(polygon)
                    else:
                        # No intersection - keep original polygon
                        new_polygons.append(polygon)

                # Update current polygons for next iteration
                current_polygons = new_polygons

            logger.info(
                f"Country {idx} split into {len(current_polygons)} regions using {split_count} boundaries"
            )

            # Create region entries from final polygons
            for i, polygon in enumerate(current_polygons):
                if polygon.is_valid and not polygon.is_empty:
                    # Create region data
                    region_data = country.copy()
                    region_data["geometry"] = polygon
                    region_data["region_id"] = f"region_{region_id:03d}"
                    region_data["original_country_id"] = idx
                    region_data["sub_region_id"] = i
                    region_data["area_km2"] = polygon.area / 1000000  # Convert to km²
                    region_data["num_boundaries_used"] = split_count

                    regions.append(region_data)
                    region_id += 1

        except Exception as e:
            logger.error(f"Error processing country {idx}: {e}")
            # Fallback: keep original geometry
            region_data = country.copy()
            region_data["region_id"] = f"region_{region_id:03d}"
            region_data["original_country_id"] = idx
            region_data["sub_region_id"] = 0
            region_data["area_km2"] = country.geometry.area / 1000000
            region_data["num_boundaries_used"] = 0
            regions.append(region_data)
            region_id += 1

    # Create GeoDataFrame from regions
    if regions:
        regions_gdf = gpd.GeoDataFrame(regions, crs=country_shapes.crs)
        logger.info(
            f"Successfully created {len(regions_gdf)} regions from {len(country_shapes)} original shapes"
        )
        return regions_gdf
    else:
        logger.error("No regions were created!")
        return country_shapes.copy()


def clean_regions(regions_gdf, min_area_threshold=1000):
    """
    Clean up regions by removing very small polygons and fixing invalid geometries.

    Args:
        regions_gdf (geopandas.GeoDataFrame): Regions to clean
        min_area_threshold (float): Minimum area threshold for keeping regions

    Returns:
        geopandas.GeoDataFrame: Cleaned regions
    """
    logger.info("Cleaning regions...")

    initial_count = len(regions_gdf)

    # Fix invalid geometries
    regions_gdf["geometry"] = regions_gdf["geometry"].buffer(0)

    # Calculate areas
    regions_gdf["area"] = regions_gdf.geometry.area

    # Remove very small regions
    regions_gdf = regions_gdf[regions_gdf["area"] > min_area_threshold]

    # Reset index
    regions_gdf = regions_gdf.reset_index(drop=True)

    logger.info(
        f"Cleaned regions: {initial_count} -> {len(regions_gdf)} (removed {initial_count - len(regions_gdf)} small regions)"
    )

    return regions_gdf


def filter_regions_with_powerplants(regions_gdf, powerplants_path):
    """
    Filter regions to keep only those that contain powerplants.
    
    Args:
        regions_gdf (geopandas.GeoDataFrame): Regions to filter
        powerplants_path (str): Path to powerplants CSV file
        
    Returns:
        geopandas.GeoDataFrame: Filtered regions containing powerplants
    """
    logger.info(f"Loading powerplants from: {powerplants_path}")
    
    try:
        # Load powerplants data
        powerplants_df = pd.read_csv(powerplants_path)
        logger.info(f"Loaded {len(powerplants_df)} powerplants")
        
        # Check for required columns
        if 'lat' not in powerplants_df.columns or 'lon' not in powerplants_df.columns:
            logger.error("Powerplants file must contain 'lat' and 'lon' columns")
            return regions_gdf
        
        # Remove rows with missing coordinates
        powerplants_df = powerplants_df.dropna(subset=['lat', 'lon'])
        logger.info(f"After removing missing coordinates: {len(powerplants_df)} powerplants")
        
        if len(powerplants_df) == 0:
            logger.warning("No powerplants with valid coordinates found")
            return regions_gdf
        
        # Create GeoDataFrame from powerplants
        powerplant_geometries = [Point(lon, lat) for lon, lat in zip(powerplants_df['lon'], powerplants_df['lat'])]
        powerplants_gdf = gpd.GeoDataFrame(
            powerplants_df, 
            geometry=powerplant_geometries, 
            crs='EPSG:4326'  # Assuming lat/lon coordinates are in WGS84
        )
        
        # Convert powerplants to same CRS as regions
        if regions_gdf.crs != powerplants_gdf.crs:
            powerplants_gdf = powerplants_gdf.to_crs(regions_gdf.crs)
            logger.info(f"Converted powerplants CRS from EPSG:4326 to {regions_gdf.crs}")
        
        # Find regions that contain powerplants
        regions_with_powerplants = []
        powerplant_count_per_region = []
        
        for idx, region in regions_gdf.iterrows():
            # Check which powerplants fall within this region
            powerplants_in_region = powerplants_gdf[powerplants_gdf.geometry.within(region.geometry)]
            
            if len(powerplants_in_region) > 0:
                regions_with_powerplants.append(region)
                powerplant_count_per_region.append(len(powerplants_in_region))
                logger.info(f"Region {region.get('region_id', idx)} contains {len(powerplants_in_region)} powerplants")
            else:
                logger.debug(f"Region {region.get('region_id', idx)} contains no powerplants - removing")
        
        if regions_with_powerplants:
            # Create new GeoDataFrame with only regions containing powerplants
            filtered_regions = gpd.GeoDataFrame(regions_with_powerplants, crs=regions_gdf.crs)
            filtered_regions['powerplant_count'] = powerplant_count_per_region
            filtered_regions = filtered_regions.reset_index(drop=True)
            
            logger.info(f"Filtered regions: {len(regions_gdf)} -> {len(filtered_regions)} (kept only regions with powerplants)")
            logger.info(f"Total powerplants in kept regions: {sum(powerplant_count_per_region)}")
            
            return filtered_regions
        else:
            logger.warning("No regions contain powerplants! Keeping all regions.")
            return regions_gdf
            
    except Exception as e:
        logger.error(f"Error filtering regions with powerplants: {e}")
        logger.warning("Continuing with unfiltered regions")
        return regions_gdf


def save_regions(regions_gdf, output_path):
    """
    Save regions to GeoJSON file.

    Args:
        regions_gdf (geopandas.GeoDataFrame): Regions to save
        output_path (str): Output file path
    """
    try:
        logger.info(f"Saving regions to: {output_path}")

        # Ensure output directory exists
        output_dir = Path(output_path).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created directory: {output_dir}")

        # Convert back to WGS84 for better map compatibility
        regions_wgs84 = regions_gdf.to_crs("EPSG:4326")

        # Save to GeoJSON in WGS84 format
        regions_wgs84.to_file(output_path, driver="GeoJSON")

        logger.info(f"Successfully saved {len(regions_gdf)} regions to {output_path}")

    except Exception as e:
        logger.error(f"Error saving regions: {e}")
        raise


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("create_region_shapes", configfiles="config/config.GB.yaml")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # load country shapes
    country_shapes = load_country_shapes(snakemake.input.country_shapes)

    # Filter out GB shapes
    country_shapes = country_shapes[country_shapes.name == "GB"]

    # load ETYS boudary lines
    boundary_lines = load_boundary_lines(snakemake.input.etys_boundary_lines)

    # Ensure boundary lines exist
    if boundary_lines.empty:
        raise ValueError("No boundary lines found in the provided ETYS data!")
    
    # create regions from boundaries
    regions = create_regions_from_boundaries(country_shapes, boundary_lines)
    logger.info(f"Created {len(regions)} initial regions")
    if len(regions) > 1:
        logger.info(
            f"- Area range: {regions.geometry.area.min() / 1000000:.1f} - {regions.geometry.area.max() / 1000000:.1f} km²"
        )
        logger.info(f"- Average area: {regions.geometry.area.mean() / 1000000:.1f} km²")

    # Clean regions with appropriate threshold (1 km²)
    min_area = 1000000  # 1 km² in square meters
    logger.info(f"\nCleaning regions (removing regions < {min_area/1000000:.0f} km²)...")
    cleaned_regions = clean_regions(regions, min_area_threshold=min_area)

    # save regions to output file
    save_regions(cleaned_regions, snakemake.output.raw_region_shapes)

    # log final summary
    logger.info("REGION CREATION SUMMARY")
    logger.info(f"ETYS boundary lines: {len(boundary_lines)}")
    logger.info(f"Regions after splitting: {len(regions)}")
