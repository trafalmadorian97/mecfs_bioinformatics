from pathlib import Path


def build_bcftools_index_command(pth: Path) -> list[str]:
    return [
        "pixi",
        "r",
        "bcftools",
        "index",
        str(pth)
    ]
