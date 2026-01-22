"""
Mapping from chromosome names found in dbsnp vcf file to standard chromosome names

"""
# List of chroosome names in dbsnp vcf file::
#
# ##contig=<ID=NC_000001.10>
# ##contig=<ID=NC_000002.11>
# ##contig=<ID=NC_000003.11>
# ##contig=<ID=NC_000004.11>
# ##contig=<ID=NC_000005.9>
# ##contig=<ID=NC_000006.11>
# ##contig=<ID=NC_000007.13>
# ##contig=<ID=NC_000008.10>
# ##contig=<ID=NC_000009.11>
# ##contig=<ID=NC_000010.10>
# ##contig=<ID=NC_000011.9>
# ##contig=<ID=NC_000012.11>
# ##contig=<ID=NC_000013.10>
# ##contig=<ID=NC_000014.8>
# ##contig=<ID=NC_000015.9>
# ##contig=<ID=NC_000016.9>
# ##contig=<ID=NC_000017.10>
# ##contig=<ID=NC_000018.9>
# ##contig=<ID=NC_000019.9>
# ##contig=<ID=NC_000020.10>
# ##contig=<ID=NC_000021.8>
# ##contig=<ID=NC_000022.10>
# ##contig=<ID=NC_000023.10>
# ##contig=<ID=NC_000024.9>
# ##contig=<ID=NC_012920.1>

import pandas as pd

NC_CHROM_MAPPING = {
    "NC_000001.10": 1,
    "NC_000002.11": 2,
    "NC_000003.10": 3,
    "NC_000004.11": 4,
    "NC_000005.9": 5,
    "NC_000006.11": 6,
    "NC_000007.13": 7,
    "NC_000008.10": 8,
    "NC_000009.11": 9,
    "NC_000010.10": 10,
    "NC_000011.9": 11,
    "NC_000012.11": 12,
    "NC_000013.10": 13,
    "NC_000014.8": 14,
    "NC_000015.9": 15,
    "NC_000016.9": 16,
    "NC_000017.10": 17,
    "NC_000018.9": 18,
    "NC_000019.9": 19,
    "NC_000020.10": 20,
    "NC_000021.8": 21,
    "NC_000022.10": 22,
    "NC_000023.10": 23,
    "NC_000024.9": 24,
    "NC_012920.1": 25,
}

NC_CHROM_TABLE = pd.DataFrame(
    {
        "nc_chrom": NC_CHROM_MAPPING.keys(),
        "int_chrom": NC_CHROM_MAPPING.values(),
    }
)
