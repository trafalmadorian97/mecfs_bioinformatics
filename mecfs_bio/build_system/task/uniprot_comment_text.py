"""
Cleaning of the free-text comment fields of a bulk UniProt export, such as "Function [CC]" and
"Subcellular location [CC]".

A comment field starts with a label (e.g. "FUNCTION:"), and its statements carry evidence code
blocks such as {ECO:0000269|PubMed:123, ECO:0000250|UniProtKB:Q8R4Y0}, where each ECO (Evidence
and Conclusion Ontology) code is optionally followed by "|" and its source. Consecutive comments
of one field are joined by "; " after a full stop; inside a comment, "; " separates items (for
example "Cell membrane; Single-pass type I membrane protein").
"""

import polars as pl

# Separator between consecutive comments of one field.
_COMMENT_SEPARATOR_PATTERN = r"\.;\s+"

# One evidence block, with the space before it.
EVIDENCE_BLOCK_PATTERN = r"\s*\{ECO:[^}]*\}"
# A block between two full stops ("Binds X. {ECO:...}."), dropped with one of them.
_EVIDENCE_SENTENCE_PATTERN = r"\.\s*\{ECO:[^}]*\}\."


def null_if_blank(expr: pl.Expr) -> pl.Expr:
    """The string with surrounding whitespace stripped, or null when nothing is left."""
    stripped = expr.str.strip_chars()
    return pl.when(stripped == "").then(None).otherwise(stripped)


def uniprot_comment_text(column: str, label: str) -> pl.Expr:
    """A UniProt comment field with its label (e.g. "FUNCTION") and evidence blocks dropped."""
    return null_if_blank(
        pl.col(column)
        .cast(pl.String)
        .str.replace_all(f"{label}:", "", literal=True)
        .str.replace_all(_EVIDENCE_SENTENCE_PATTERN, ".")
        .str.replace_all(EVIDENCE_BLOCK_PATTERN, "")
        .str.replace_all(_COMMENT_SEPARATOR_PATTERN, ". ")
        .str.replace_all(r"\s+", " ")
    )
