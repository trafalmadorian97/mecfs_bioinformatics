#!/usr/bin/env bash
# Does filesystem access matter for unit-test time?
# 1. Profile one full serial run with cProfile (output: suite.prof, analysed by profile_categories.py).
# Serial timings in fs_cost.log were lost to a since-fixed bc dependency; see gc_and_fs_benchmark.sh.
# 2. Time full runs with pytest temp dirs on the ext4 disk vs on tmpfs (/dev/shm, RAM),
#    serial and xdist, interleaved, two repetitions each.
# Usage (from repo root): bash experiments/claude/test_speed/fs_cost_benchmark.sh | tee experiments/claude/test_speed/fs_cost.log
set -u
OUT=experiments/claude/test_speed
COMMON="-p no:testmon -p no:cacheprovider -q --typeguard-packages=mecfs_bio test_mecfs_bio/unit"
CAPS="OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_MAX_THREADS=2 POLARS_MAX_THREADS=2"

echo "=== cProfile serial run ==="
start=$(date +%s.%N)
pixi r python -m cProfile -o $OUT/suite.prof -m pytest $COMMON | tail -2
awk -v s=$start -v e=$(date +%s.%N) 'BEGIN{printf "profiled_wall_s=%.1f\n", e-s}'

run() {  # label, tmp root, extra pytest args, env prefix
  local label=$1 root=$2 extra=$3 envp=$4
  rm -rf "$root/fsbench" && mkdir -p "$root/fsbench/tmp"
  local s=$(date +%s.%N)
  env $envp TMPDIR=$root/fsbench/tmp pixi r python -m pytest $COMMON --basetemp=$root/fsbench/base $extra 2>&1 | tail -1
  awk -v s=$s -v e=$(date +%s.%N) -v l=$label 'BEGIN{printf "RESULT %s wall_s=%.1f\n", l, e-s}'
  rm -rf "$root/fsbench"
}
DISK=$HOME/.cache/fsbench_disk
mkdir -p $DISK
for rep in 1 2; do
  run serial_disk  $DISK     ""      ""
  run serial_tmpfs /dev/shm  ""      ""
  run xdist_disk   $DISK     "-n 8"  "$CAPS"
  run xdist_tmpfs  /dev/shm  "-n 8"  "$CAPS"
done
rmdir $DISK 2>/dev/null
