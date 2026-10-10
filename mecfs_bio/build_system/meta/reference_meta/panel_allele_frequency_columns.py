"""Which allele-frequency column of a reference panel table holds which ancestry."""

from collections.abc import Sequence

from attrs import frozen

from mecfs_bio.constants.allele_frequency_panel_constants import (
    PanelAncestry,
    panel_af_col,
)


@frozen(slots=True)
class PanelAncestryColumn:
    ancestry: PanelAncestry
    column: str


@frozen(slots=True)
class PanelAlleleFrequencyColumns:
    """The ancestries a panel carries, each with its allele-frequency column.

    Entries are a tuple so the declaration stays hashable on a frozen Task's meta.
    """

    entries: tuple[PanelAncestryColumn, ...]

    def __attrs_post_init__(self):
        ancestries = self.ancestries
        assert ancestries, "a panel must declare at least one ancestry"
        assert len(set(ancestries)) == len(ancestries), (
            f"duplicate panel ancestries: {ancestries}"
        )

    @property
    def ancestries(self) -> list[PanelAncestry]:
        return [entry.ancestry for entry in self.entries]

    def column_for(self, ancestry: PanelAncestry) -> str:
        matches = [entry.column for entry in self.entries if entry.ancestry == ancestry]
        assert matches, (
            f"the panel has no allele frequency for ancestry {ancestry!r}; "
            f"available: {self.ancestries}"
        )
        return matches[0]

    @classmethod
    def single(
        cls, ancestry: PanelAncestry, column: str
    ) -> "PanelAlleleFrequencyColumns":
        return cls((PanelAncestryColumn(ancestry=ancestry, column=column),))

    @classmethod
    def prefixed(
        cls, ancestries: Sequence[PanelAncestry]
    ) -> "PanelAlleleFrequencyColumns":
        """One AF_<ancestry> column per ancestry, as multi-ancestry panels name them."""
        return cls(
            tuple(
                PanelAncestryColumn(ancestry=ancestry, column=panel_af_col(ancestry))
                for ancestry in ancestries
            )
        )
