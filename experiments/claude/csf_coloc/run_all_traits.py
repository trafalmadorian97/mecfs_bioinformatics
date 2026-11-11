"""Run the CSF coloc of the non-AD traits one after another in one process (never two runners
at once). Usage: pixi r python experiments/claude/csf_coloc/run_all_traits.py [name...]"""

import sys
import time

from mecfs_bio.analysis import csf_coloc_analysis as a

RUNS = {
    "decode_me": a.run_decode_me_csf_coloc,
    "ra": a.run_ra_seropositive_csf_coloc,
    "ldl": a.run_willer_ldl_csf_coloc,
    "ibd": a.run_liu_et_al_2023_ibd_csf_coloc,
    "scz": a.run_pgc_2022_sch_csf_coloc,
    "pd": a.run_nalls_2019_parkinsons_csf_coloc,
    "als": a.run_van_rheenen_2021_als_csf_coloc,
    "decode_me_exploratory": a.run_decode_me_csf_coloc_exploratory,
}


def main() -> None:
    for name in sys.argv[1:] or list(RUNS):
        start = time.time()
        print(f"=== START {name}", flush=True)
        RUNS[name]()
        print(f"=== DONE {name} in {(time.time() - start) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
