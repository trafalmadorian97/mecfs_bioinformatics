#!/usr/bin/env bash
# Follow-up to fs_cost_benchmark.sh (whose serial timings were lost: bc is not installed).
# Times full serial runs with temp dirs on disk vs tmpfs, and with/without gc.freeze() after
# collection; also xdist with/without gc.freeze(). Reports total time spent in gc.collect().
# Usage (repo root): bash experiments/claude/test_speed/gc_and_fs_benchmark.sh > experiments/claude/test_speed/gc_and_fs.log 2>&1
set -u
COMMON="-p no:testmon -p no:cacheprovider -q -s -p gc_timing_plugin --typeguard-packages=mecfs_bio test_mecfs_bio/unit"
CAPS="OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_MAX_THREADS=2 POLARS_MAX_THREADS=2"
DISK=$HOME/.cache/fsbench_disk
run() {  # label, tmp root, extra pytest args, env prefix
  local label=$1 root=$2 extra=$3 envp=$4
  rm -rf "$root/fsbench" && mkdir -p "$root/fsbench/tmp"
  local s=$(date +%s.%N)
  env $envp PYTHONPATH=experiments/claude/test_speed TMPDIR=$root/fsbench/tmp \
    pixi r python -m pytest $COMMON --basetemp=$root/fsbench/base $extra > "$root/fsbench.out" 2>&1
  local e=$(date +%s.%N)
  grep -a "passed\|failed" "$root/fsbench.out" | tail -1
  grep -a "GCSTATS" "$root/fsbench.out" | awk '{for(i=1;i<=NF;i++){split($i,a,"=");if(a[1]=="seconds")s+=a[2];if(a[1]=="calls")c+=a[2]}} END{printf "gc_calls=%d gc_seconds_summed_over_workers=%.1f\n",c,s}'
  awk -v s=$s -v e=$e -v l=$label 'BEGIN{printf "RESULT %s wall_s=%.1f\n", l, e-s}'
  rm -rf "$root/fsbench" "$root/fsbench.out"
}
mkdir -p $DISK
for rep in 1 2; do
  run serial_disk         $DISK    ""     ""
  run serial_tmpfs        /dev/shm ""     ""
  run serial_disk_freeze  $DISK    ""     "GC_FREEZE=1"
  run xdist_disk          $DISK    "-n 8" "$CAPS"
  run xdist_disk_freeze   $DISK    "-n 8" "$CAPS GC_FREEZE=1"
done
rmdir $DISK 2>/dev/null
