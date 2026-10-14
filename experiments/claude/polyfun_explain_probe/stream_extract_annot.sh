#!/usr/bin/env bash
# THROWAWAY spike helper. Stream the polyfun baseline-LF 2.2.UKB annotation
# tarball from S3 and extract only *.annot.parquet members, logging every
# extracted member name. Aborts naturally when the stream ends; the caller
# kills it once enough annotation files have landed.
set -uo pipefail
URL="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.polyfun.tar.gz"
OUT_DIR="$(dirname "$0")/annot_out"
LOG="$(dirname "$0")/logs/extract.log"
mkdir -p "$OUT_DIR"
echo "[$(date -Is)] starting stream extract of *.annot.parquet" | tee -a "$LOG"
# --wildcards matches the member basename pattern anywhere in the path.
curl -s "$URL" | tar -xzv -C "$OUT_DIR" --wildcards '*.annot.parquet' 2>&1 | tee -a "$LOG"
echo "[$(date -Is)] stream ended" | tee -a "$LOG"
