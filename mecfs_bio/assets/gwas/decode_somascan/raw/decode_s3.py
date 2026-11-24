"""
Where deCODE serves the Eldjarn et al. 2023 SomaScan summary statistics: an S3-compatible
bucket at deCODE, read with keys deCODE issues on request. The keys live in the AWS
credentials profile DECODE_CREDENTIALS_PROFILE (aws configure --profile decode-download), never
in the repository.
"""

from mecfs_bio.build_system.wf.s3_access import S3Endpoint, S3ObjectRef
from mecfs_bio.constants.decode_somascan_constants import DecodeSomascanGroup

DECODE_CREDENTIALS_PROFILE = "decode-download"
DECODE_S3_ENDPOINT = S3Endpoint(
    url="https://s3-ext.decode.is",
    region="us-east-1",
    credentials_profile=DECODE_CREDENTIALS_PROFILE,
)
DECODE_BUCKET = "largescaleplasma-2023"
DECODE_PREFIX: dict[DecodeSomascanGroup, str] = {
    "smp": "final_somascan_smp/",
    "raw": "final_somascan_raw/",
}
DECODE_VARIANT_INFO_REF = S3ObjectRef(
    endpoint=DECODE_S3_ENDPOINT,
    bucket=DECODE_BUCKET,
    key="variant_info_iceland.txt.gz",
)
