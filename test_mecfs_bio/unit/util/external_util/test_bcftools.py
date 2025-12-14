from pathlib import Path

from mecfs_bio.util.external_util.bcftools import build_bcftools_index_command
from mecfs_bio.util.subproc.run_command import execute_command


def test_bcftools_indexing_command(tmp_path: Path):
    """
    Test that we can construct an index using the bcftools indexing command
    """
    vcf_file_path = tmp_path / "output.vcf"
    compressed_path = tmp_path / "output.vcf.gz"
    index_path = tmp_path / "output.vcf.gz.csi"
    #
    assert not compressed_path.exists()
    assert not index_path.exists()
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	1000	.	A	G	100	PASS	DP=10	GT	0/1
chr1	1050	rs999	C	T	.	.	DP=20	GT	1/1"""
    with open(vcf_file_path, "w") as f:
        f.write(vcf_content)
    execute_command(["pixi", "r", "bgzip", str(vcf_file_path)])
    assert compressed_path.exists()
    execute_command(
        build_bcftools_index_command(
            compressed_path,
        )
    )
    assert index_path.exists()
