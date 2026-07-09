"""
Unit tests for converting POPs' gene_annot_jun10.txt into a MAGMA gene-location file.
"""

from pathlib import Path

from mecfs_bio.build_system.task.pops.pops_magma_gene_loc_task import (
    write_magma_gene_loc,
)

_ANNOT = (
    "ENSGID\tNAME\tCHR\tSTART\tEND\tTSS\n"
    "ENSG00000186092\tOR4F5\t1\t69091\t70008\t69091\n"  # TSS == START -> +
    "ENSG00000188976\tNOC2L\t1\t879584\t894689\t894689\n"  # TSS == END -> -
)


def test_write_magma_gene_loc(tmp_path: Path):
    """Produces MAGMA-format rows (ENSGID CHR START END STRAND NAME) with strand
    inferred from the TSS."""
    annot = tmp_path / "gene_annot_jun10.txt"
    annot.write_text(_ANNOT, encoding="utf-8")
    out = tmp_path / "out.gene.loc"

    num_genes = write_magma_gene_loc(annot, out)

    assert num_genes == 2
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert rows[0] == ["ENSG00000186092", "1", "69091", "70008", "+", "OR4F5"]
    assert rows[1] == ["ENSG00000188976", "1", "879584", "894689", "-", "NOC2L"]


def test_write_magma_gene_loc_rejects_bad_header(tmp_path: Path):
    """A source table with an unexpected header fails fast."""
    annot = tmp_path / "bad.txt"
    annot.write_text("GENE\tCHR\tSTART\n", encoding="utf-8")
    out = tmp_path / "out.gene.loc"

    try:
        write_magma_gene_loc(annot, out)
    except AssertionError as err:
        assert "header" in str(err)
    else:
        raise AssertionError("expected AssertionError for bad header")
