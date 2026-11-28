"""Download of DAS files from an S3-compatible bucket."""
from __future__ import annotations

# --- Standard library --------------------------------------------------------
import datetime
import logging
import os

logger = logging.getLogger(__name__)


def get_data_from_S3(
    bucket_name: str, endpoint_url: str, starttime: datetime.datetime,
    endtime: datetime.datetime, local_path: str, unsigned: bool = False,
) -> bool:
    """
    Download DAS raw files from an S3-compatible endpoint for a single day.

    Parameters:
    - bucket_name: name of the S3 bucket.
    - endpoint_url: S3 endpoint URL.
    - starttime, endtime: datetime.datetime objects defining the interval (must be same day).
    - local_path: local directory where files will be saved (organized by hour).
    - unsigned: if True, use unsigned (anonymous) S3 access.

    Returns:
    - True on success, False on failure or invalid input.
    """
    # Optional dependency: imported lazily so the module works without boto3
    import boto3
    from botocore.config import Config
    from botocore import UNSIGNED

    # 1. Multi-day Guard: Warn and exit if start and end are on different days
    if starttime.date() != endtime.date():
        logger.warning(f"WARNING: starttime ({starttime.date()}) and endtime ({endtime.date()}) "
              "are on different days. This function only processes single-day requests.")
        return False

    # 2. Setup S3 Resource
    if unsigned:
        s3_res = boto3.resource('s3', endpoint_url=endpoint_url, config=Config(signature_version=UNSIGNED))
    else:
        s3_res = boto3.resource('s3', endpoint_url=endpoint_url)

    bucket = s3_res.Bucket(bucket_name)

    try:
        dir_structure = starttime.strftime('%Y%m%d/dphi/')
        h5start = starttime.strftime('%Y%m%d/dphi/%H%M%S.h5')
        h5end = endtime.strftime('%Y%m%d/dphi/%H%M%S.h5')
        
        files_to_download = []
        last_seen_key = None
        
        # Check if we should skip the buffer file (Midnight check)
        is_midnight = (starttime.hour == 0 and starttime.minute == 0 and starttime.second == 0)

        for obj in bucket.objects.filter(Prefix=dir_structure):
            if obj.key > h5end:
                break
            
            if h5start <= obj.key <= h5end:
                # Add buffer file ONLY if not midnight AND it's the first file found
                if not files_to_download and last_seen_key and not is_midnight:
                    files_to_download.append(last_seen_key)
                
                files_to_download.append(obj.key)
            
            last_seen_key = obj.key

        # 3. Execution
        if not files_to_download:
            logger.warning("No files found for the specified time range.")
            return True

        target_dir = os.path.join(local_path, dir_structure)
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        for key in files_to_download:
            target = os.path.join(target_dir, os.path.basename(key))
            if os.path.exists(target):
                logger.info(f"File {target} already exists in {target_dir}. Skipping download.")
                continue
            logger.info(f"Downloading: {key} to {target}")
            bucket.download_file(key, target)

        logger.info(f"Successfully downloaded {len(files_to_download)} files.")
        return True

    except Exception as e:
        logger.error(f"Error: {e}")
        return False
