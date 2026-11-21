#!/usr/bin/env bash
# End-to-end test: download N whole SomaScan files concurrently with aria2 (16 connections each,
# via presigned URLs), report aggregate throughput, and verify each against the md5 deCODE
# stores in the object's metadata (x-amz-meta-md5chksum, base64). Needs AWS_* creds in env.
# Usage: aria2_whole_files.sh <n_files> <dest_dir>
set -euo pipefail
N=$1; DEST=$2; EP=https://s3-ext.decode.is; B=largescaleplasma-2023
mkdir -p "$DEST"
mapfile -t KEYS < <(aws --endpoint-url $EP s3 ls s3://$B/final_somascan_smp/ | awk '{print "final_somascan_smp/"$4}' | sed -n "20,$((19+N))p")
: > "$DEST/input.txt"
for k in "${KEYS[@]}"; do
  echo "$(aws --endpoint-url $EP s3 presign s3://$B/$k --expires-in 7200)" >> "$DEST/input.txt"
  echo "  out=$(basename "$k")" >> "$DEST/input.txt"
done
t0=$(date +%s.%N)
pixi r -e download-env aria2c -q -i "$DEST/input.txt" -d "$DEST" -j "$N" -x 16 -s 16 -k 8M \
  --file-allocation=none --max-tries=8 --retry-wait=5 --timeout=60 --allow-overwrite=true
t1=$(date +%s.%N)
bytes=$(du -cb "$DEST"/*.txt.gz | tail -1 | cut -f1)
python3 -c "print(f'{$N} files {$bytes/1e9:.2f} GB in {$t1-$t0:.0f}s = {$bytes/1e6/($t1-$t0):.1f} MB/s aggregate')"
for k in "${KEYS[@]}"; do
  want=$(aws --endpoint-url $EP s3api head-object --bucket $B --key "$k" --query Metadata.md5chksum --output text | base64 -d | xxd -p)
  got=$(md5sum "$DEST/$(basename "$k")" | cut -d' ' -f1)
  [[ "$want" == "$got" ]] && echo "md5 OK  $(basename "$k")" || echo "md5 BAD $(basename "$k") want=$want got=$got"
done
