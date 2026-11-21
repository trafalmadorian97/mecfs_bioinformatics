#!/usr/bin/env bash
# Measure aggregate deCODE S3 throughput vs number of concurrent HTTP Range connections.
# Each connection fetches a distinct 16 MB range via a presigned URL with curl (no AWS CLI in
# the data path). Requires AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY in the environment.
# Usage: parallel_scaling.sh "1 4 8 16 32"
set -euo pipefail
EP=https://s3-ext.decode.is
B=largescaleplasma-2023
KEYS=(final_somascan_smp/Proteomics_SMP_PC0_10000_28_CRYBB2_CRBB2_10032022.txt.gz
      final_somascan_smp/Proteomics_SMP_PC0_10001_7_RAF1_c_Raf_10032022.txt.gz
      final_somascan_raw/Proteomics_PC0_10003_15_ZNF41_ZNF41_07082019.txt.gz
      final_somascan_raw/Proteomics_PC0_10006_25_ELK1_ELK1_07082019.txt.gz)
URLS=()
for k in "${KEYS[@]}"; do URLS+=("$(aws --endpoint-url $EP s3 presign s3://$B/$k --expires-in 7200)"); done
CHUNK=$((16*1024*1024))
for n in ${1:-1 4 8 16 32}; do
  start=$(date +%s.%N)
  for i in $(seq 0 $((n-1))); do
    u=${URLS[$((i % ${#URLS[@]}))]}
    off=$(( ((37*i) % 900) * 1024 * 1024 + RANDOM*4096 ))
    curl -sS -o /dev/null -r $off-$((off+CHUNK-1)) -w "%{size_download}\n" "$u" &
  done > "${TMPDIR:-/tmp}/claude_scaling_sizes_$n.txt"
  wait
  end=$(date +%s.%N)
  total=$(awk '{s+=$1} END {print s}' "${TMPDIR:-/tmp}/claude_scaling_sizes_$n.txt")
  python3 -c "print(f'n={$n:3d} bytes={$total/1e6:.0f}MB seconds={$end-$start:.1f} aggregate={$total/1e6/($end-$start):.1f}MB/s per_conn={$total/1e6/($end-$start)/$n:.2f}MB/s')"
  rm -f "${TMPDIR:-/tmp}/claude_scaling_sizes_$n.txt"
done
