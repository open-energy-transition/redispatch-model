# SPDX-FileCopyrightText:  gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT

"""
Retrieve Generation Unit Unavailability Data from ENTSO-E Transparency Platform

This script retrieves generation unit unavailability data from the ENTSO-E API
for the UK bidding zones and processes it for use in PyPSA modeling.

API Documentation: 15.1.A&B Unavailability of Generation Units
"""

import io
import logging
import os
import time
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Optional

import pandas as pd
import requests
from dotenv import load_dotenv

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


class ENTSOEUnavailabilityAPI:
    """
    Interface to ENTSO-E Transparency Platform API for unavailability data
    """

    def __init__(self, api_key: str):
        """
        Initialize the API client

        Args:
            api_key: ENTSO-E API key for authentication
        """
        self.api_key = api_key
        self.base_url = "https://web-api.tp.entsoe.eu/api"
        self.session = requests.Session()
        self.rate_limit_delay = 2.0  # Seconds between requests

    def _make_request(self, params: dict[str, Any]) -> requests.Response:
        """
        Make API request with rate limiting and error handling

        Args:
            params: Query parameters for the API request

        Returns:
            Response object
        """
        # Add API key to parameters
        params["securityToken"] = self.api_key

        logger.debug(f"Making API request with params: {params}")

        try:
            # Rate limiting
            time.sleep(self.rate_limit_delay)

            response = self.session.get(self.base_url, params=params, timeout=30)
            logger.debug(f"API response status: {response.status_code}")
            logger.debug(f"API response length: {len(response.content)} bytes")

            if response.status_code != 200:
                logger.error(f"API returned status {response.status_code}")

            response.raise_for_status()
            return response

        except requests.exceptions.RequestException as e:
            logger.error(f"API request failed: {e}")
            raise

    def get_unavailability_data(
        self,
        bidding_zone: str,
        period_start: datetime,
        period_end: datetime,
        business_type: Optional[str] = None,
        doc_status: Optional[str] = None,
        registered_resource: Optional[str] = None,
        offset: int = 0,
    ) -> bytes:
        """
        Retrieve unavailability data from ENTSO-E API

        Args:
            bidding_zone: EIC code for bidding zone (e.g., '10YGB----------A' for GB)
            period_start: Start datetime for the query period
            period_end: End datetime for the query period
            business_type: Optional business type ('A53' for planned, 'A54' for forced)
            doc_status: Optional document status ('A05', 'A09', 'A13')
            registered_resource: Optional EIC code of specific generation unit
            offset: Offset for pagination (0-4800)

        Returns:
            Raw response content as bytes (could be XML or ZIP)
        """
        params = {
            "documentType": "A80",  # Generation unavailability
            "BiddingZone_Domain": bidding_zone,
            "periodStart": period_start.strftime("%Y%m%d%H%M"),
            "periodEnd": period_end.strftime("%Y%m%d%H%M"),
        }

        # Add optional parameters
        if business_type:
            params["BusinessType"] = business_type
        if doc_status:
            params["DocStatus"] = doc_status
        if registered_resource:
            params["RegisteredResource"] = registered_resource
        if offset > 0:
            params["offset"] = offset

        response = self._make_request(params)
        return response.content


def generate_weekly_periods(
    start_date: datetime, end_date: datetime, max_days: int = 7
) -> list[tuple]:
    """
    Generate weekly periods for API requests to stay within API limits

    Args:
        start_date: Overall start date
        end_date: Overall end date
        max_days: Maximum days per request (default: 7)

    Returns:
        List of (start, end) datetime tuples for each period
    """
    periods = []
    current_start = start_date

    while current_start < end_date:
        current_end = min(current_start + timedelta(days=max_days), end_date)
        periods.append((current_start, current_end))
        current_start = current_end

    logger.info(
        f"Split {start_date.date()} to {end_date.date()} into {len(periods)} weekly periods"
    )
    return periods


def retrieve_data_by_periods(
    api_client,
    bidding_zone: str,
    periods: list[tuple],
    business_type: str,
    doc_status: str,
    xml_save_dir: str,
) -> list[pd.DataFrame]:
    """
    Retrieve data for multiple periods and return list of DataFrames

    Args:
        api_client: ENTSOEUnavailabilityAPI instance
        bidding_zone: EIC code for bidding zone
        periods: List of (start, end) datetime tuples
        business_type: Business type code
        doc_status: Document status code
        xml_save_dir: Directory to save XML files

    Returns:
        List of DataFrames, one per period
    """
    all_dataframes = []
    successful_periods = 0

    # Log overall progress
    logger.info(
        f"Retrieving {len(periods)} periods from {periods[0][0].date()} to {periods[-1][1].date()}"
    )

    for i, (period_start, period_end) in enumerate(periods, 1):
        # Show progress every 10 periods or at start/end
        if i == 1 or i == len(periods) or i % 10 == 0:
            logger.info(
                f"Progress: {i}/{len(periods)} periods ({(i / len(periods) * 100):.0f}%)"
            )

        try:
            # Add period identifier to XML save directory
            period_xml_dir = f"{xml_save_dir}/period_{period_start.strftime('%Y%m%d')}_{period_end.strftime('%Y%m%d')}"

            response_content = api_client.get_unavailability_data(
                bidding_zone=bidding_zone,
                period_start=period_start,
                period_end=period_end,
                business_type=business_type,
                doc_status=doc_status,
            )

            # Parse data and save XML files
            df = parse_unavailability_zip(response_content, period_xml_dir)

            if not df.empty:
                # Add period metadata
                df["request_period_start"] = period_start
                df["request_period_end"] = period_end
                df["period_batch"] = i
                all_dataframes.append(df)
                successful_periods += 1

            # Rate limiting between requests
            if i < len(periods):
                time.sleep(2)  # 2 second delay between requests

        except Exception as e:
            # Only log actual errors, not "no data" cases
            if "File is not a zip file" not in str(e):
                logger.warning(
                    f"Period {i} ({period_start.date()} to {period_end.date()}): {e}"
                )
            continue

    logger.info(
        f"Retrieved data from {successful_periods}/{len(periods)} periods ({(successful_periods / len(periods) * 100):.0f}% success rate)"
    )
    return all_dataframes


def merge_period_dataframes(dataframes: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Merge DataFrames from multiple periods and remove duplicates

    Args:
        dataframes: List of DataFrames to merge

    Returns:
        Combined DataFrame with duplicates removed
    """
    if not dataframes:
        logger.warning("No dataframes to merge")
        return pd.DataFrame()

    # Combine all dataframes
    combined_df = pd.concat(dataframes, ignore_index=True)
    logger.info(
        f"Combined {len(dataframes)} periods into {len(combined_df)} total records"
    )

    # Remove duplicates based on key fields
    # Use resource_mrid and start_time/end_time as the primary key for deduplication
    # This focuses on actual outage periods rather than data reporting periods
    key_columns = ["resource_mrid", "start_time", "end_time"]

    # Only use columns that exist in the dataframe
    existing_key_columns = [col for col in key_columns if col in combined_df.columns]

    # If start_time/end_time not available, fall back to resource_mrid and business_type
    if (
        "start_time" not in existing_key_columns
        and "resource_mrid" in combined_df.columns
    ):
        key_columns = ["resource_mrid", "business_type"]
        existing_key_columns = [
            col for col in key_columns if col in combined_df.columns
        ]

    if existing_key_columns:
        initial_count = len(combined_df)
        combined_df = combined_df.drop_duplicates(
            subset=existing_key_columns, keep="first"
        )
        removed_count = initial_count - len(combined_df)

        if removed_count > 0:
            logger.info(
                f"Removed {removed_count} duplicate records, {len(combined_df)} records remaining"
            )

    # Sort by start time and resource name
    if "start_time" in combined_df.columns:
        combined_df = combined_df.sort_values(["start_time", "resource_name"])
    else:
        combined_df = combined_df.sort_values(["resource_name"])

    return combined_df


def parse_unavailability_zip(
    zip_content: bytes, save_xml_dir: Optional[str] = None
) -> pd.DataFrame:
    """
    Parse ZIP file containing XML unavailability data

    Args:
        zip_content: ZIP file content as bytes
        save_xml_dir: Optional directory to save extracted XML files

    Returns:
        DataFrame with parsed unavailability data
    """
    outages = []

    try:
        with zipfile.ZipFile(io.BytesIO(zip_content)) as zf:
            xml_files = [f for f in zf.namelist() if f.endswith(".xml")]
            # Only log if we actually find files
            if xml_files:
                logger.debug(f"Processing ZIP file with {len(xml_files)} XML files")

            # Create save directory if specified
            if save_xml_dir:
                save_path = Path(save_xml_dir)
                save_path.mkdir(parents=True, exist_ok=True)
                # Only log if we have files to save
                if xml_files:
                    logger.debug(f"Saving XML files to: {save_path}")

            for filename in xml_files:
                try:
                    xml_content = zf.read(filename).decode("utf-8")

                    # Save XML file if directory specified
                    if save_xml_dir:
                        xml_file_path = save_path / filename
                        with open(xml_file_path, "w", encoding="utf-8") as f:
                            f.write(xml_content)
                        # Removed individual file save logging for brevity

                    file_outages = parse_xml_content(xml_content, filename)
                    outages.extend(file_outages)

                except Exception as e:
                    logger.error(f"Error processing {filename}: {e}")
                    continue

    except zipfile.BadZipFile as e:
        logger.error(f"Invalid ZIP file: {e}")
        return pd.DataFrame()

    if not outages:
        # Don't log "no data" as it's common
        return pd.DataFrame()

    df = pd.DataFrame(outages)
    if len(df) > 0:
        logger.debug(f"Parsed {len(df)} outage records from ZIP file")
    return df


def parse_xml_content(xml_content: str, filename: str = "") -> list[dict]:
    """
    Parse individual XML content for unavailability data

    Args:
        xml_content: XML content as string
        filename: Optional filename for logging

    Returns:
        List of outage dictionaries
    """
    outages = []

    try:
        root = ET.fromstring(xml_content)
        namespaces = {"ns": "urn:iec62325.351:tc57wg16:451-6:outagedocument:3:0"}

        # Find all TimeSeries elements
        timeseries_elements = root.findall(".//ns:TimeSeries", namespaces)

        for ts in timeseries_elements:
            outage_data = {}

            # TimeSeries information
            outage_data["business_type"] = _get_text(
                ts, ".//ns:businessType", namespaces
            )
            outage_data["bidding_zone"] = _get_text(
                ts, ".//ns:biddingZone_Domain.mRID", namespaces
            )

            # Resource information
            outage_data["resource_mrid"] = _get_text(
                ts, ".//ns:production_RegisteredResource.mRID", namespaces
            )
            outage_data["resource_name"] = _get_text(
                ts, ".//ns:production_RegisteredResource.name", namespaces
            )
            outage_data["resource_location"] = _get_text(
                ts, ".//ns:production_RegisteredResource.location.name", namespaces
            )
            outage_data["resource_type"] = _get_text(
                ts, ".//ns:production_RegisteredResource.pSRType.psrType", namespaces
            )

            # Power system resource details
            outage_data["psr_mrid"] = _get_text(
                ts,
                ".//ns:production_RegisteredResource.pSRType.powerSystemResources.mRID",
                namespaces,
            )
            outage_data["psr_name"] = _get_text(
                ts,
                ".//ns:production_RegisteredResource.pSRType.powerSystemResources.name",
                namespaces,
            )
            nominal_power = _get_text(
                ts,
                ".//ns:production_RegisteredResource.pSRType.powerSystemResources.nominalP",
                namespaces,
            )
            outage_data["nominal_power_mw"] = (
                float(nominal_power) if nominal_power else None
            )

            # Time periods from TimeSeries level
            start_date = _get_text(ts, ".//ns:start_DateAndOrTime.date", namespaces)
            start_time = _get_text(ts, ".//ns:start_DateAndOrTime.time", namespaces)
            end_date = _get_text(ts, ".//ns:end_DateAndOrTime.date", namespaces)
            end_time = _get_text(ts, ".//ns:end_DateAndOrTime.time", namespaces)

            if start_date and start_time:
                outage_data["start_time"] = pd.to_datetime(f"{start_date} {start_time}")
            if end_date and end_time:
                outage_data["end_time"] = pd.to_datetime(f"{end_date} {end_time}")

            # Process Available_Period elements
            available_periods = ts.findall(".//ns:Available_Period", namespaces)

            if available_periods:
                for period in available_periods:
                    period_data = outage_data.copy()

                    # Extract data points
                    points = period.findall(".//ns:Point", namespaces)
                    if points:
                        quantities = []
                        for point in points:
                            qty = _get_text(point, ".//ns:quantity", namespaces)
                            if qty:
                                quantities.append(float(qty))

                        if quantities:
                            period_data["min_available_mw"] = min(quantities)
                            period_data["max_available_mw"] = max(quantities)
                            period_data["avg_available_mw"] = sum(quantities) / len(
                                quantities
                            )
                            period_data["num_data_points"] = len(quantities)

                            # Calculate unavailable capacity
                            if outage_data["nominal_power_mw"]:
                                period_data["min_unavailable_mw"] = max(
                                    0, outage_data["nominal_power_mw"] - max(quantities)
                                )
                                period_data["max_unavailable_mw"] = max(
                                    0, outage_data["nominal_power_mw"] - min(quantities)
                                )

                    outages.append(period_data)
            else:
                # No Available_Period elements, just add the TimeSeries data
                outages.append(outage_data)

    except ET.ParseError as e:
        logger.error(f"Failed to parse XML {filename}: {e}")
    except Exception as e:
        logger.error(f"Error processing XML {filename}: {e}")

    return outages


def _get_text(element, xpath: str, namespaces: dict[str, str]) -> Optional[str]:
    """
    Helper function to safely extract text from XML element

    Args:
        element: XML element
        xpath: XPath expression
        namespaces: XML namespaces

    Returns:
        Text content or None
    """
    found = element.find(xpath, namespaces)
    return found.text if found is not None else None


def process_unavailability_data(
    df: pd.DataFrame, output_file: str, bidding_zone: str, config: dict
) -> None:
    """
    Process and save unavailability data

    Args:
        df: DataFrame with raw unavailability data
        output_file: Path to save processed data
        bidding_zone: Bidding zone identifier
        config: Snakemake configuration dictionary
    """
    if df.empty:
        logger.warning("No data to process")
        # Create empty CSV with headers
        empty_df = pd.DataFrame(
            columns=[
                "resource_mrid",
                "resource_name",
                "resource_location",
                "resource_type",
                "business_type",
                "start_time",
                "end_time",
                "nominal_power_mw",
                "min_available_mw",
                "max_available_mw",
                "avg_available_mw",
                "min_unavailable_mw",
                "max_unavailable_mw",
                "bidding_zone",
                "business_type_desc",
            ]
        )
        empty_df.to_csv(output_file, index=False)
        return

    # Add bidding zone information
    df["bidding_zone"] = bidding_zone

    # Clean and process data
    entso_e_config = config.get("unavailability", {}).get("entso_e", {})
    business_type_map = {
        entso_e_config.get("business_types", {}).get(
            "planned", "A53"
        ): "Planned maintenance",
        entso_e_config.get("business_types", {}).get(
            "forced", "A54"
        ): "Forced unavailability",
    }
    df["business_type_desc"] = df["business_type"].map(business_type_map)

    # Sort by start time
    df = df.sort_values(["start_time", "resource_name"])

    # Select only the columns we want to keep in the output
    output_columns = [
        "resource_mrid",
        "resource_name",
        "resource_location",
        "resource_type",
        "business_type",
        "start_time",
        "end_time",
        "nominal_power_mw",
        "min_available_mw",
        "max_available_mw",
        "avg_available_mw",
        "min_unavailable_mw",
        "max_unavailable_mw",
        "bidding_zone",
        "business_type_desc",
    ]

    # Filter to only include columns that exist in the dataframe
    columns_to_keep = [col for col in output_columns if col in df.columns]
    df_filtered = df[columns_to_keep]

    # Save processed data
    df_filtered.to_csv(output_file, index=False)
    logger.info(f"Saved {len(df_filtered)} records to {output_file}")

    # Print summary statistics
    logger.info("Data summary:")
    if (
        "start_time" in df_filtered.columns
        and not df_filtered["start_time"].isnull().all()
    ):
        logger.info(
            f"  Date range: {df_filtered['start_time'].min()} to {df_filtered['end_time'].max()}"
        )
    logger.info(f"  Total outage records: {len(df_filtered)}")
    logger.info(f"  Unique generation units: {df_filtered['resource_name'].nunique()}")

    if (
        "nominal_power_mw" in df_filtered.columns
        and not df_filtered["nominal_power_mw"].isnull().all()
    ):
        logger.info(
            f"  Total installed capacity: {df_filtered['nominal_power_mw'].sum():.1f} MW"
        )

    # Business type breakdown
    if "business_type_desc" in df_filtered.columns:
        business_breakdown = (
            df_filtered.groupby("business_type_desc")
            .agg({"resource_name": "nunique", "resource_mrid": "count"})
            .round(1)
        )
        business_breakdown.columns = ["Unique_Units", "Total_Records"]
        logger.info(f"  Breakdown by type:\\n{business_breakdown}")


def test_api_connection(api_client, zone_code: str) -> bool:
    """
    Test API connection with a small date range

    Args:
        api_client: ENTSOEUnavailabilityAPI instance
        zone_code: EIC code for bidding zone

    Returns:
        True if API connection works
    """
    try:
        # Test with a recent date range
        test_start = datetime.now() - timedelta(days=30)
        test_end = datetime.now() - timedelta(days=29)

        logger.info(
            f"Testing API connection with date range: {test_start.date()} to {test_end.date()}"
        )

        response_content = api_client.get_unavailability_data(
            bidding_zone=zone_code,
            period_start=test_start,
            period_end=test_end,
            business_type="A53",  # Planned maintenance
        )

        logger.info(f"Test API response length: {len(response_content)} bytes")

        # Check if it's a ZIP file
        is_zip = response_content.startswith(b"PK")
        logger.info(f"Response is ZIP format: {is_zip}")

        return True

    except Exception as e:
        logger.error(f"API connection test failed: {e}")
        return False


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_unavailability_data")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Load environment variables from .env file
    load_dotenv()

    # Get API key from environment variable
    api_key = os.getenv("ENTSO_E_API_KEY")

    if not api_key:
        raise ValueError(
            "ENTSO-E API key not found. Please set ENTSO_E_API_KEY in your .env file or environment variables.\\n"
            "You can get an API key from: https://transparency.entsoe.eu/usrm/user/createPublicApiUser.do"
        )

    # Get date parameters from config or use defaults
    unavailability_config = snakemake.config.get("unavailability", {})
    start_date = pd.to_datetime(unavailability_config.get("start_date", "2023-01-01"))
    end_date = pd.to_datetime(unavailability_config.get("end_date", "2023-12-31"))

    bidding_zones = unavailability_config.get("bidding_zones", ["GB"])
    business_types = unavailability_config.get("business_types", ["planned", "forced"])
    max_request_days = unavailability_config.get("max_request_days", 7)

    # Get mapping dictionaries from config
    entso_e_config = unavailability_config.get("entso_e", {})
    uk_bidding_zones = entso_e_config.get("bidding_zones", {})
    business_types_map = entso_e_config.get("business_types", {})
    doc_status_map = entso_e_config.get("doc_status", {})

    # Initialize API client
    api_client = ENTSOEUnavailabilityAPI(api_key)

    # Test API connection first
    logger.info("Testing API connection...")
    if not test_api_connection(api_client, uk_bidding_zones.get("GB")):
        raise ConnectionError(
            "API connection test failed. Check your API key and network connection."
        )

    # Process each bidding zone
    for zone in bidding_zones:
        if zone not in uk_bidding_zones:
            logger.warning(f"Unknown bidding zone: {zone}")
            continue

        zone_code = uk_bidding_zones[zone]
        logger.info(f"Retrieving unavailability data for {zone} ({zone_code})")

        # Process each business type
        for business_type in business_types:
            if business_type not in business_types_map:
                logger.warning(f"Unknown business type: {business_type}")
                continue

            business_code = business_types_map[business_type]
            logger.info(f"Processing {business_type} outages ({business_code})")

            # Split the date range into weekly periods
            periods = generate_weekly_periods(
                start_date, end_date, max_days=max_request_days
            )

            # Create XML storage directory
            xml_save_dir = f"data/gb-model/outage_data/{zone.lower()}_{business_type}"

            # Retrieve data for all periods
            all_dataframes = retrieve_data_by_periods(
                api_client=api_client,
                bidding_zone=zone_code,
                periods=periods,
                business_type=business_code,
                doc_status=doc_status_map["active"],
                xml_save_dir=xml_save_dir,
            )

            # Merge all periods into a single DataFrame
            df = merge_period_dataframes(all_dataframes)
            logger.info(f"Final merged DataFrame has {len(df)} rows")

            # Generate output filename
            output_file = snakemake.output[f"{business_type}_unavailability"]

            # Process and save data
            process_unavailability_data(df, output_file, zone, snakemake.config)

            # Log XML storage information across all periods
            xml_base_dir = Path(xml_save_dir)
            if xml_base_dir.exists():
                # Count XML files across all period subdirectories
                total_xml_count = len(list(xml_base_dir.rglob("*.xml")))
                total_size = sum(f.stat().st_size for f in xml_base_dir.rglob("*.xml"))

                period_dirs = [
                    d
                    for d in xml_base_dir.iterdir()
                    if d.is_dir() and d.name.startswith("period_")
                ]

                logger.info(
                    f"Stored {total_xml_count} XML files ({total_size / 1024 / 1024:.1f} MB) across {len(period_dirs)} periods in {xml_base_dir}"
                )

    logger.info("Unavailability data retrieval completed")
