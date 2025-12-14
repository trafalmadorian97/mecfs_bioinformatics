"""
Functions for using the bcftools utiltiies.
"""

from pathlib import Path


def build_bcftools_index_command(pth: Path) -> list[str]:
    """
    Returns a command line command to construct bcftools index on a given compressed vcf file.
    """
    return ["pixi", "r", "bcftools", "index", str(pth)]
