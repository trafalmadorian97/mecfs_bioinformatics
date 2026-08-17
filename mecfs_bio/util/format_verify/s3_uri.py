"""
Pragmatic validation of s3:// URIs (bucket plus optional key).

Catches realistic mistakes at graph-construction time -- a missing s3:// scheme,
embedded whitespace, or a bucket name that violates the core S3 naming rules --
rather than implementing every corner of the bucket-naming specification.
Accepts s3://mecfs-bio-remote-exec-scratch/remote-exec-scratch and s3://bucket
(key optional); rejects a bare path with no scheme, s3://Bad_Bucket, and "".
"""

import re

_SCHEME = "s3://"
# Charset + endpoint rule: lowercase alnum / dot / hyphen, starting and ending
# with an alphanumeric. Length is checked separately.
_BUCKET_CHARSET_RE = re.compile(r"[a-z0-9](?:[a-z0-9.-]*[a-z0-9])?")
# A bucket name must not be formatted as an IPv4 address.
_IPV4_RE = re.compile(r"[0-9]{1,3}(?:\.[0-9]{1,3}){3}")


def is_valid_s3_bucket_name(name: str) -> bool:
    """Return whether name satisfies the core S3 bucket-naming rules.

    Enforces the length bound (3-63), the lowercase alnum/./- charset with
    alphanumeric endpoints, no consecutive dots, and no IPv4-address form.
    Some rarer restrictions (reserved prefixes/suffixes such as xn-- or
    -s3alias) are intentionally not checked -- this is a pragmatic guard, not
    the full specification.
    """
    if not 3 <= len(name) <= 63:
        return False
    if not _BUCKET_CHARSET_RE.fullmatch(name):
        return False
    if ".." in name:
        return False
    if _IPV4_RE.fullmatch(name):
        return False
    return True


def is_valid_s3_uri(uri: str) -> bool:
    """Return whether uri is a plausibly valid s3://bucket[/key] address.

    Requires the s3:// scheme, no embedded whitespace, and a bucket name that
    passes is_valid_s3_bucket_name. The key is optional and, being effectively
    unrestricted in S3, is only required to be free of whitespace (already
    guaranteed by the whole-string check).
    """
    if not uri or any(char.isspace() for char in uri):
        return False
    if not uri.startswith(_SCHEME):
        return False
    bucket, _, _key = uri[len(_SCHEME) :].partition("/")
    return is_valid_s3_bucket_name(bucket)
