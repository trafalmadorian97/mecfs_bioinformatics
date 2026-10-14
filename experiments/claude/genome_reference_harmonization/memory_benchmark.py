"""
V4: compare peak anonymous memory of genome-reference harmonization with gwaslab harmonization.

Each scenario runs in its own process. Acceptance: the all-chromosome peak is at most
1.5 times the chromosome-1-and-2 peak, i.e. peak memory follows the largest chromosome
rather than total rows.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.memory_benchmark \
    2>&1 | tee experiments/claude/genome_reference_harmonization/memory_benchmark.log
"""

import json

from mecfs_bio.util.subproc.run_command import execute_command

SCENARIOS = ["genome_reference_chr1_2", "genome_reference_all", "gwaslab"]
ACCEPTABLE_RATIO = 1.5


def main() -> None:
    peaks: dict[str, float] = {}
    for scenario in SCENARIOS:
        output = execute_command(
            [
                "python",
                "-m",
                "experiments.claude.genome_reference_harmonization.measure_harmonization_memory",
                scenario,
            ]
        )
        result = json.loads(output.strip().splitlines()[-1])
        print(result)
        peaks[scenario] = result["peak_rss_anon_gib"]
    ratio = peaks["genome_reference_all"] / peaks["genome_reference_chr1_2"]
    print(f"all / chr1-2 peak ratio: {ratio:.2f} (acceptable <= {ACCEPTABLE_RATIO})")
    print(
        f"gwaslab / genome-reference peak ratio: {peaks['gwaslab'] / peaks['genome_reference_all']:.2f}"
    )
    assert ratio <= ACCEPTABLE_RATIO, (
        "peak memory grows with total rows; investigate before switching"
    )


if __name__ == "__main__":
    main()
