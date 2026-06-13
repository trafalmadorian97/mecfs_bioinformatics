"""
One-off: replace the hand-pasted ``??? info "Variant List..."`` admonition
blocks in the DecodeME SUSIE fine-mapping docs with ``markdown_table`` macro
calls pointing at the programmatically generated credible-set tables.

Each block is: a line ``??? info "<title>"`` followed by a blank line and a
4-space-indented markdown table.  We remove the whole block (header line +
following blank/indented lines) and substitute a single macro call.  The macro
re-creates the collapsible automatically because these tables are large.
"""

from pathlib import Path

DOCS = Path("docs/Analysis/ME_CFS/DecodeME/Fine_Mapping/SUSIE")

# file -> list of (admonition title, .mdx asset_id)
REPLACEMENTS = {
    "a_Chr1_173M_174M_Locus.md": [
        (
            "Variant List",
            "decode_mechr1_173500000_174500000_palindromes_keep_susie_base_convert_cs_to_markdown",
        ),
    ],
    "c_Chr6_97M_99M_Locus.md": [
        (
            "Variant List (L=10)",
            "decode_mechr6_97500000_99000000_palindromes_keep_susie_base_convert_cs_to_markdown",
        ),
        (
            "Variant List (L=2)",
            "decode_mechr6_97500000_99000000_palindromes_keep_susie_2_convert_cs_to_markdown",
        ),
    ],
    "d_Chr15_54M_55M_Locus.md": [
        (
            "Variant List (L=10)",
            "decode_mechr15_54500000_55500000_palindromes_keep_susie_base_convert_cs_to_markdown",
        ),
    ],
    "e_Chr17_50M_51M_Locus.md": [
        (
            "Variant List (L=10)",
            "decode_mechr17_50000000_51000000_palindromes_keep_susie_base_convert_cs_to_markdown",
        ),
    ],
    "f_Chr20_47M_48M_Locus.md": [
        (
            "Variant List (L=10)",
            "decode_mechr20_47000000_48200000_palindromes_keep_susie_base_convert_cs_to_markdown",
        ),
    ],
}


def macro_call(title: str, asset_id: str) -> str:
    return f'{{{{ markdown_table("docs/_figs/{asset_id}.mdx", title="{title}") }}}}'


def replace_block(lines: list[str], title: str, asset_id: str) -> list[str]:
    header = f'??? info "{title}"'
    for i, line in enumerate(lines):
        if line.rstrip("\n") == header:
            # Consume the header and every following blank or 4-space-indented
            # line (the admonition body).
            j = i + 1
            while j < len(lines):
                stripped = lines[j].rstrip("\n")
                if stripped == "" or stripped.startswith("    "):
                    j += 1
                else:
                    break
            new = lines[:i] + [macro_call(title, asset_id) + "\n", "\n"] + lines[j:]
            return new
    raise ValueError(f"admonition header not found: {header!r}")


def main():
    for filename, repls in REPLACEMENTS.items():
        path = DOCS / filename
        lines = path.read_text().splitlines(keepends=True)
        for title, asset_id in repls:
            lines = replace_block(lines, title, asset_id)
        path.write_text("".join(lines))
        print(f"updated {path} ({len(repls)} table(s))")


if __name__ == "__main__":
    main()
