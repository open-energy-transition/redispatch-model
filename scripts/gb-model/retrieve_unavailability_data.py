# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Retrieve Generation Unit Unavailability Data from ENTSO-E Transparency Platform

This script retrieves generation unit unavailability data from the ENTSO-E API
for the UK bidding zones and processes it for use in PyPSA modeling.

API Documentation: 15.1.A&B Unavailability of Generation Units
"""

import requests
import pandas as pd
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import logging
from pathlib import Path
from typing import Optional, Dict, List, Any
import time
import os
from dotenv import load_dotenv
import zipfile
import io
import shutil

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

# UK bidding zone codes
UK_BIDDING_ZONES = {
    "GB": "10YGB----------A",  # Great Britain
}

# Business types for outage classification
BUSINESS_TYPES = {
    "planned": "A53",      # Planned maintenance
    "forced": "A54"        # Forced unavailability (unplanned outage)
}

# Document status codes
DOC_STATUS = {
    "active": "A05",       # Active
    "cancelled": "A09",    # Cancelled
    "withdrawn": "A13"     # Withdrawn
}


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
        self.rate_limit_delay = 1.0  # Seconds between requests
        
    def _make_request(self, params: Dict[str, Any]) -> requests.Response:
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
        offset: int = 0
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


def generate_weekly_periods(start_date: datetime, end_date: datetime, max_days: int = 7) -> List[tuple]:
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
    
    logger.info(f"Split {start_date.date()} to {end_date.date()} into {len(periods)} weekly periods")
    return periods


def retrieve_data_by_periods(
    api_client,
    bidding_zone: str,
    periods: List[tuple],
    business_type: str,
    doc_status: str,
    xml_save_dir: str
) -> List[pd.DataFrame]:
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
    
    for i, (period_start, period_end) in enumerate(periods, 1):
        logger.info(f"Retrieving period {i}/{len(periods)}: {period_start.date()} to {period_end.date()}")
        
        try:
            # Add period identifier to XML save directory
            period_xml_dir = f"{xml_save_dir}/period_{period_start.strftime('%Y%m%d')}_{period_end.strftime('%Y%m%d')}"
            
            response_content = api_client.get_unavailability_data(
                bidding_zone=bidding_zone,
                period_start=period_start,
                period_end=period_end,
                business_type=business_type,
                doc_status=doc_status
            )
            
            logger.debug(f"Retrieved {len(response_content)} bytes for period {i}")
            
            # Parse data and save XML files
            df = parse_unavailability_zip(response_content, period_xml_dir)
            
            if not df.empty:
                # Add period metadata
                df['request_period_start'] = period_start
                df['request_period_end'] = period_end
                df['period_batch'] = i
                all_dataframes.append(df)
                logger.info(f"Period {i}: Parsed {len(df)} records")
            else:
                logger.info(f"Period {i}: No data found")
            
            # Rate limiting between requests
            if i < len(periods):
                time.sleep(2)  # 2 second delay between requests
                
        except Exception as e:
            logger.error(f"Failed to retrieve period {i} ({period_start.date()} to {period_end.date()}): {e}")
            continue
    
    logger.info(f"Successfully retrieved {len(all_dataframes)} periods out of {len(periods)}")
    return all_dataframes


def merge_period_dataframes(dataframes: List[pd.DataFrame]) -> pd.DataFrame:
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
    logger.info(f"Combined {len(dataframes)} periods into {len(combined_df)} total records")
    
    # Remove duplicates based on key fields
    # Use document_mrid and start_time/end_time as the primary key for deduplication
    # This focuses on actual outage periods rather than data reporting periods
    key_columns = ['document_mrid', 'timeseries_mrid', 'start_time', 'end_time']
    
    # Only use columns that exist in the dataframe
    existing_key_columns = [col for col in key_columns if col in combined_df.columns]
    
    # If start_time/end_time not available, fall back to period_start/period_end
    if 'start_time' not in existing_key_columns and 'period_start' in combined_df.columns:
        key_columns = ['document_mrid', 'timeseries_mrid', 'period_start', 'period_end']
        existing_key_columns = [col for col in key_columns if col in combined_df.columns]
    
    if existing_key_columns:
        initial_count = len(combined_df)
        combined_df = combined_df.drop_duplicates(subset=existing_key_columns, keep='first')
        removed_count = initial_count - len(combined_df)
        
        if removed_count > 0:
            logger.info(f"Removed {removed_count} duplicate records, {len(combined_df)} records remaining")
    
    # Sort by start time and resource name
    if 'start_time' in combined_df.columns:
        combined_df = combined_df.sort_values(['start_time', 'resource_name'])
    elif 'period_start' in combined_df.columns:
        combined_df = combined_df.sort_values(['period_start', 'resource_name'])
    
    return combined_df


def parse_unavailability_zip(zip_content: bytes, save_xml_dir: Optional[str] = None) -> pd.DataFrame:
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
            xml_files = [f for f in zf.namelist() if f.endswith('.xml')]
            logger.info(f"Processing ZIP file with {len(xml_files)} XML files")
            
            # Create save directory if specified
            if save_xml_dir:
                save_path = Path(save_xml_dir)
                save_path.mkdir(parents=True, exist_ok=True)
                logger.info(f"Saving XML files to: {save_path}")
            
            for filename in xml_files:
                try:
                    xml_content = zf.read(filename).decode('utf-8')
                    
                    # Save XML file if directory specified
                    if save_xml_dir:
                        xml_file_path = save_path / filename
                        with open(xml_file_path, 'w', encoding='utf-8') as f:
                            f.write(xml_content)
                        logger.debug(f"Saved {filename} to {xml_file_path}")
                    
                    file_outages = parse_xml_content(xml_content, filename)
                    outages.extend(file_outages)
                    
                except Exception as e:
                    logger.error(f"Error processing {filename}: {e}")
                    continue
                    
    except zipfile.BadZipFile as e:
        logger.error(f"Invalid ZIP file: {e}")
        return pd.DataFrame()
    
    if not outages:
        logger.warning("No outage data found in ZIP file")
        return pd.DataFrame()
    
    df = pd.DataFrame(outages)
    logger.info(f"Parsed {len(df)} outage records from ZIP file")
    return df


def parse_xml_content(xml_content: str, filename: str = "") -> List[Dict]:
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
        namespaces = {'ns': 'urn:iec62325.351:tc57wg16:451-6:outagedocument:3:0'}
        
        # Extract document-level information
        doc_mrid = _get_text(root, './/ns:mRID', namespaces)
        doc_revision = _get_text(root, './/ns:revisionNumber', namespaces)
        doc_status = _get_text(root, './/ns:docStatus/ns:value', namespaces)
        
        # Find all TimeSeries elements
        timeseries_elements = root.findall('.//ns:TimeSeries', namespaces)
        
        for ts in timeseries_elements:
            outage_data = {}
            
            # Document information
            outage_data['document_mrid'] = doc_mrid
            outage_data['revision_number'] = doc_revision
            outage_data['doc_status'] = doc_status
            outage_data['source_file'] = filename
            
            # TimeSeries information
            outage_data['timeseries_mrid'] = _get_text(ts, './/ns:mRID', namespaces)
            outage_data['business_type'] = _get_text(ts, './/ns:businessType', namespaces)
            outage_data['bidding_zone'] = _get_text(ts, './/ns:biddingZone_Domain.mRID', namespaces)
            
            # Resource information
            outage_data['resource_mrid'] = _get_text(ts, './/ns:production_RegisteredResource.mRID', namespaces)
            outage_data['resource_name'] = _get_text(ts, './/ns:production_RegisteredResource.name', namespaces)
            outage_data['resource_location'] = _get_text(ts, './/ns:production_RegisteredResource.location.name', namespaces)
            outage_data['resource_type'] = _get_text(ts, './/ns:production_RegisteredResource.pSRType.psrType', namespaces)
            
            # Power system resource details
            outage_data['psr_mrid'] = _get_text(ts, './/ns:production_RegisteredResource.pSRType.powerSystemResources.mRID', namespaces)
            outage_data['psr_name'] = _get_text(ts, './/ns:production_RegisteredResource.pSRType.powerSystemResources.name', namespaces)
            nominal_power = _get_text(ts, './/ns:production_RegisteredResource.pSRType.powerSystemResources.nominalP', namespaces)
            outage_data['nominal_power_mw'] = float(nominal_power) if nominal_power else None
            
            # Time periods from TimeSeries level
            start_date = _get_text(ts, './/ns:start_DateAndOrTime.date', namespaces)
            start_time = _get_text(ts, './/ns:start_DateAndOrTime.time', namespaces)
            end_date = _get_text(ts, './/ns:end_DateAndOrTime.date', namespaces)
            end_time = _get_text(ts, './/ns:end_DateAndOrTime.time', namespaces)
            
            if start_date and start_time:
                outage_data['start_time'] = pd.to_datetime(f"{start_date} {start_time}")
            if end_date and end_time:
                outage_data['end_time'] = pd.to_datetime(f"{end_date} {end_time}")
            
            # Process Available_Period elements
            available_periods = ts.findall('.//ns:Available_Period', namespaces)
            
            if available_periods:
                for period in available_periods:
                    period_data = outage_data.copy()
                    
                    # Period-specific time interval
                    interval_start = _get_text(period, './/ns:timeInterval/ns:start', namespaces)
                    interval_end = _get_text(period, './/ns:timeInterval/ns:end', namespaces)
                    
                    if interval_start:
                        period_data['period_start'] = pd.to_datetime(interval_start)
                    if interval_end:
                        period_data['period_end'] = pd.to_datetime(interval_end)
                    
                    # Resolution
                    period_data['resolution'] = _get_text(period, './/ns:resolution', namespaces)
                    
                    # Extract data points
                    points = period.findall('.//ns:Point', namespaces)
                    if points:
                        quantities = []
                        for point in points:
                            qty = _get_text(point, './/ns:quantity', namespaces)
                            if qty:
                                quantities.append(float(qty))
                        
                        if quantities:
                            period_data['min_available_mw'] = min(quantities)
                            period_data['max_available_mw'] = max(quantities)
                            period_data['avg_available_mw'] = sum(quantities) / len(quantities)
                            period_data['num_data_points'] = len(quantities)
                            
                            # Calculate unavailable capacity
                            if outage_data['nominal_power_mw']:
                                period_data['min_unavailable_mw'] = max(0, outage_data['nominal_power_mw'] - max(quantities))
                                period_data['max_unavailable_mw'] = max(0, outage_data['nominal_power_mw'] - min(quantities))
                    
                    outages.append(period_data)
            else:
                # No Available_Period elements, just add the TimeSeries data
                outages.append(outage_data)
                
    except ET.ParseError as e:
        logger.error(f"Failed to parse XML {filename}: {e}")
    except Exception as e:
        logger.error(f"Error processing XML {filename}: {e}")
    
    return outages


def _get_text(element, xpath: str, namespaces: Dict[str, str]) -> Optional[str]:
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
    df: pd.DataFrame,
    output_file: str,
    bidding_zone: str
) -> None:
    """
    Process and save unavailability data
    
    Args:
        df: DataFrame with raw unavailability data
        output_file: Path to save processed data
        bidding_zone: Bidding zone identifier
    """
    if df.empty:
        logger.warning("No data to process")
        # Create empty CSV with headers
        empty_df = pd.DataFrame(columns=[
            'document_mrid', 'resource_mrid', 'resource_name', 'resource_location',
            'resource_type', 'business_type', 'doc_status',
            'start_time', 'end_time', 'period_start', 'period_end',
            'nominal_power_mw', 'min_available_mw', 'max_available_mw', 'avg_available_mw',
            'min_unavailable_mw', 'max_unavailable_mw',
            'bidding_zone'
        ])
        empty_df.to_csv(output_file, index=False)
        return
    
    # Add bidding zone information
    df['bidding_zone'] = bidding_zone
    
    # Clean and process data
    df['business_type_desc'] = df['business_type'].map({
        'A53': 'Planned maintenance',
        'A54': 'Forced unavailability'
    })
    
    df['doc_status_desc'] = df['doc_status'].map({
        'A05': 'Active',
        'A09': 'Cancelled',
        'A13': 'Withdrawn'
    })
    
    # Sort by start time
    df = df.sort_values(['start_time', 'resource_name'])
    
    # Save processed data
    df.to_csv(output_file, index=False)
    logger.info(f"Saved {len(df)} records to {output_file}")
    
    # Print summary statistics
    logger.info("Data summary:")
    if 'start_time' in df.columns and not df['start_time'].isnull().all():
        logger.info(f"  Date range: {df['start_time'].min()} to {df['end_time'].max()}")
    logger.info(f"  Total outage records: {len(df)}")
    logger.info(f"  Unique generation units: {df['resource_name'].nunique()}")
    
    if 'nominal_power_mw' in df.columns and not df['nominal_power_mw'].isnull().all():
        logger.info(f"  Total installed capacity: {df['nominal_power_mw'].sum():.1f} MW")
    
    # Business type breakdown
    if 'business_type_desc' in df.columns:
        business_breakdown = df.groupby('business_type_desc').agg({
            'resource_name': 'nunique',
            'document_mrid': 'count'
        }).round(1)
        business_breakdown.columns = ['Unique_Units', 'Total_Records']
        logger.info(f"  Breakdown by type:\\n{business_breakdown}")


def create_xml_index(xml_dir: str) -> pd.DataFrame:
    """
    Create an index of stored XML files with metadata
    
    Args:
        xml_dir: Directory containing XML files
        
    Returns:
        DataFrame with XML file metadata
    """
    xml_path = Path(xml_dir)
    if not xml_path.exists():
        return pd.DataFrame()
    
    xml_files = list(xml_path.glob('*.xml'))
    if not xml_files:
        return pd.DataFrame()
    
    index_data = []
    
    for xml_file in xml_files:
        try:
            # Extract metadata from filename
            filename = xml_file.name
            
            # Parse basic file info
            file_info = {
                'filename': filename,
                'filepath': str(xml_file),
                'size_bytes': xml_file.stat().st_size,
                'modified_time': datetime.fromtimestamp(xml_file.stat().st_mtime)
            }
            
            # Try to extract period from filename
            if 'PLANNED_UNAVAIL_OF_GENERATION_UNITS_' in filename:
                parts = filename.split('PLANNED_UNAVAIL_OF_GENERATION_UNITS_')[1].replace('.xml', '')
                if '-' in parts:
                    start_str, end_str = parts.split('-')
                    try:
                        file_info['period_start'] = pd.to_datetime(start_str, format='%Y%m%d%H%M')
                        file_info['period_end'] = pd.to_datetime(end_str, format='%Y%m%d%H%M')
                    except:
                        pass
            
            index_data.append(file_info)
            
        except Exception as e:
            logger.warning(f"Error indexing {xml_file}: {e}")
            continue
    
    if index_data:
        df = pd.DataFrame(index_data)
        # Save index file
        index_file = xml_path / 'xml_index.csv'
        df.to_csv(index_file, index=False)
        logger.info(f"Created XML index with {len(df)} files: {index_file}")
        return df
    
    return pd.DataFrame()


def cleanup_xml_storage(base_dir: str = "data/gb-model/outage_data", max_age_days: int = 30):
    """
    Clean up old XML files to manage storage space
    
    Args:
        base_dir: Base directory for XML storage
        max_age_days: Maximum age of files to keep
    """
    base_path = Path(base_dir)
    if not base_path.exists():
        return
    
    cutoff_time = datetime.now() - timedelta(days=max_age_days)
    deleted_count = 0
    deleted_size = 0
    
    for xml_file in base_path.rglob('*.xml'):
        try:
            if datetime.fromtimestamp(xml_file.stat().st_mtime) < cutoff_time:
                file_size = xml_file.stat().st_size
                xml_file.unlink()
                deleted_count += 1
                deleted_size += file_size
        except Exception as e:
            logger.warning(f"Error deleting {xml_file}: {e}")
    
    if deleted_count > 0:
        logger.info(f"Cleaned up {deleted_count} old XML files ({deleted_size/1024/1024:.1f} MB)")


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
        
        logger.info(f"Testing API connection with date range: {test_start.date()} to {test_end.date()}")
        
        response_content = api_client.get_unavailability_data(
            bidding_zone=zone_code,
            period_start=test_start,
            period_end=test_end,
            business_type="A53",  # Planned maintenance
        )
        
        logger.info(f"Test API response length: {len(response_content)} bytes")
        
        # Check if it's a ZIP file
        is_zip = response_content.startswith(b'PK')
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
    start_date = pd.to_datetime(snakemake.config.get("unavailability_start_date", "2023-01-01"))
    end_date = pd.to_datetime(snakemake.config.get("unavailability_end_date", "2023-12-31"))
    
    bidding_zones = snakemake.config.get("unavailability_bidding_zones", ["GB"])
    business_types = snakemake.config.get("unavailability_business_types", ["planned", "forced"])
    max_request_days = snakemake.config.get("max_request_days", 7)
    
    logger.info(f"Configuration loaded:")
    logger.info(f"  Date range: {start_date.date()} to {end_date.date()}")
    logger.info(f"  Bidding zones: {bidding_zones}")
    logger.info(f"  Business types: {business_types}")
    logger.info(f"  Max request days: {max_request_days}")
    
    # Initialize API client
    api_client = ENTSOEUnavailabilityAPI(api_key)
    
    # Test API connection first
    logger.info("Testing API connection...")
    if not test_api_connection(api_client, UK_BIDDING_ZONES["GB"]):
        raise ConnectionError("API connection test failed. Check your API key and network connection.")
    
    # Process each bidding zone
    for zone in bidding_zones:
        if zone not in UK_BIDDING_ZONES:
            logger.warning(f"Unknown bidding zone: {zone}")
            continue
            
        zone_code = UK_BIDDING_ZONES[zone]
        logger.info(f"Retrieving unavailability data for {zone} ({zone_code})")
        
        # Process each business type
        for business_type in business_types:
            if business_type not in BUSINESS_TYPES:
                logger.warning(f"Unknown business type: {business_type}")
                continue
                
            business_code = BUSINESS_TYPES[business_type]
            logger.info(f"Processing {business_type} outages ({business_code})")
            
            try:
                # Split the date range into weekly periods
                periods = generate_weekly_periods(start_date, end_date, max_days=max_request_days)
                
                # Create XML storage directory
                xml_save_dir = f"data/gb-model/outage_data/{zone.lower()}_{business_type}"
                
                # Retrieve data for all periods
                all_dataframes = retrieve_data_by_periods(
                    api_client=api_client,
                    bidding_zone=zone_code,
                    periods=periods,
                    business_type=business_code,
                    doc_status=DOC_STATUS["active"],
                    xml_save_dir=xml_save_dir
                )
                
                # Merge all periods into a single DataFrame
                df = merge_period_dataframes(all_dataframes)
                logger.info(f"Final merged DataFrame has {len(df)} rows")
                
                # Generate output filename
                output_file = snakemake.output[f"{zone.lower()}_{business_type}_unavailability"]
                
                # Process and save data
                process_unavailability_data(df, output_file, zone)
                
                # Log XML storage information across all periods
                xml_base_dir = Path(xml_save_dir)
                if xml_base_dir.exists():
                    # Count XML files across all period subdirectories
                    total_xml_count = len(list(xml_base_dir.rglob('*.xml')))
                    total_size = sum(f.stat().st_size for f in xml_base_dir.rglob('*.xml'))
                    
                    period_dirs = [d for d in xml_base_dir.iterdir() if d.is_dir() and d.name.startswith('period_')]
                    
                    logger.info(f"Stored {total_xml_count} XML files ({total_size/1024/1024:.1f} MB) across {len(period_dirs)} periods in {xml_base_dir}")
                    
                    # Create XML index for the overall collection
                    create_xml_index(xml_save_dir)
                
            except Exception as e:
                logger.error(f"Failed to process {zone} {business_type} data: {e}")
                logger.exception("Full traceback:")
                continue
    
    # Clean up old XML files if configured
    cleanup_days = snakemake.config.get("xml_cleanup_days", 30)
    if cleanup_days > 0:
        cleanup_xml_storage(max_age_days=cleanup_days)
    
    logger.info("Unavailability data retrieval completed")