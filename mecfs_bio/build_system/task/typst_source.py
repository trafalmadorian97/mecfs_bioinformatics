"""
A Typst diagram source: exactly one of a committed .typ file path or an inline
string of Typst code.
"""

from pathlib import Path

from attrs import frozen

_SOURCE_FILENAME = "diagram.typ"


@frozen
class TypstSource:
    """
    Holds exactly one of a path to a .typ file or a literal string of Typst
    code. Enforcing exactly-one at construction makes the invalid states (both
    set, neither set, missing file) unrepresentable.
    """

    path: Path | None = None
    code: str | None = None

    def __attrs_post_init__(self):
        assert (self.path is None) != (
            self.code is None
        ), "TypstSource requires exactly one of path or code"
        if self.path is not None:
            assert (
                self.path.is_file()
            ), f"Typst source file does not exist: {self.path}"

    def resolve_to_path(self, scratch_dir: Path) -> Path:
        """
        Return a concrete .typ path for the source. For the path variant this is
        the existing file; for the code variant the string is written into
        scratch_dir and that new path is returned.
        """
        if self.path is not None:
            return self.path
        assert self.code is not None  # guaranteed by __attrs_post_init__
        target = scratch_dir / _SOURCE_FILENAME
        target.write_text(self.code)
        return target
