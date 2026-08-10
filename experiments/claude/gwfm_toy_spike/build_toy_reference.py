"""One-shot generator for the committed GWFM known-truth toy reference (spike artifact).

Run manually and offline; its OUTPUTS are committed under
test_mecfs_bio/system/test_data/gwfm_toy so the system test needs no Docker-based
generation at test time. This script:

  1. runs the proven synthesis gen_toy_gwfm.py into a temp dir (toy.bed/bim/fam,
     blocks.txt, toy.ma, annot.txt, genemap.txt);
  2. builds the gctb:test image from docker/gctb and runs --make-block-ldm inside it to
     produce the blockwise LD reference folder ldm13M/ (no local gctb binary assumed);
  3. packages the reference with the EXACT inner names GctbFineMapTask expects
     (ukbEUR_13M_FullLDM.zip -> ldm13M/, annot_baseline2.2_13M.zip -> a single
     annot_baseline2.2_13M.txt, gene_map_hg38_hg19.txt);
  4. writes the three reference files into test_data/gwfm_toy/reference/ and toy.ma one
     level up into test_data/gwfm_toy/.

See experiments/claude/design_specs/2026-08-09-gwfm-toy-recipe.md for the proven recipe
and the crash root-causes that shaped the synthetic data.

Usage: pixi r python experiments/claude/gwfm_toy_spike/build_toy_reference.py
"""

import os
import shutil
import tempfile
import zipfile
from pathlib import Path

from mecfs_bio.build_system.task.sbayesrc import gctb_gwfm_constants as c
from mecfs_bio.util.subproc.run_command import execute_command

_HERE = Path(__file__).resolve().parent
_REPO_ROOT = _HERE.parents[2]
_GEN_SCRIPT = _HERE / "gen_toy_gwfm.py"
_DOCKER_CONTEXT = _REPO_ROOT / "docker" / "gctb"
_IMAGE_TAG = "gctb:test"
_TOY_DATA = _REPO_ROOT / "test_mecfs_bio" / "system" / "test_data" / "gwfm_toy"
_REFERENCE_DIR = _TOY_DATA / "reference"


def _zip_dir(root: Path, dir_name: str, out_zip: Path) -> None:
    """Zip root/dir_name so the archive unzips back to dir_name/ (single-root archive)."""
    src = root / dir_name
    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(p for p in src.rglob("*") if p.is_file()):
            zf.write(path, arcname=str(path.relative_to(root)))


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)

        execute_command(["pixi", "r", "python", str(_GEN_SCRIPT), str(work)])

        execute_command(
            [
                "docker",
                "build",
                "--build-arg",
                f"GCTB_URL={c.GCTB_BINARY_URL}",
                "--build-arg",
                f"GCTB_SHA256={c.GCTB_BINARY_SHA256}",
                "-t",
                _IMAGE_TAG,
                str(_DOCKER_CONTEXT),
            ]
        )

        # Build the blockwise LD reference inside the gctb image (recon: >=2 blocks writes
        # the folder-level snp.info/ldm.info/rsq0.5.pwld implicitly; --thread 1 for
        # determinism).
        execute_command(
            [
                "docker",
                "run",
                "--rm",
                # Run as the host user so the container-written ldm13M/ is owned by us
                # (otherwise TemporaryDirectory cleanup fails on root-owned files).
                "-u",
                f"{os.getuid()}:{os.getgid()}",
                "-v",
                f"{work}:/w",
                "-w",
                "/w",
                _IMAGE_TAG,
                "gctb",
                "--bfile",
                "toy",
                "--make-block-ldm",
                "--block-info",
                "blocks.txt",
                "--out",
                c.GWFM_LDM_DIR_NAME,
                "--thread",
                "1",
            ]
        )

        _REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

        _zip_dir(work, c.GWFM_LDM_DIR_NAME, _REFERENCE_DIR / c.GWFM_LDM_ZIP_NAME)

        annot_txt = work / c.GWFM_ANNOT_FILE_NAME
        shutil.copyfile(work / "annot.txt", annot_txt)
        with zipfile.ZipFile(
            _REFERENCE_DIR / c.GWFM_ANNOT_ZIP_NAME, "w", zipfile.ZIP_DEFLATED
        ) as zf:
            zf.write(annot_txt, arcname=c.GWFM_ANNOT_FILE_NAME)

        shutil.copyfile(work / "genemap.txt", _REFERENCE_DIR / c.GWFM_GENE_MAP_FILE_NAME)
        shutil.copyfile(work / "toy.ma", _TOY_DATA / "toy.ma")

    print(f"wrote reference to {_REFERENCE_DIR}")
    for path in sorted(_REFERENCE_DIR.iterdir()):
        print(f"  {path.name}: {path.stat().st_size} bytes")
    print(f"wrote toy.ma: {(_TOY_DATA / 'toy.ma').stat().st_size} bytes")


if __name__ == "__main__":
    main()
