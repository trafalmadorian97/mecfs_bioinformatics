#!/usr/bin/env bash
# Download ONE SomaScan file with many parallel range requests, using s5cmd and rclone (neither
# has aria2's 16-connections-per-server cap), and verify the md5 deCODE stores in the object
# metadata. Needs AWS_* creds in env. Run via:
#   pixi exec --spec s5cmd --spec rclone -- single_file_many_connections.sh <dest_dir> "<conns>"
set -euo pipefail
DEST=$1; CONNS=${2:-"32 64 128"}; EP=https://s3-ext.decode.is; B=largescaleplasma-2023
K=final_somascan_smp/Proteomics_SMP_PC0_10076_1_AP4M1_AP4M1_10032022.txt.gz
mkdir -p "$DEST"
want=$(aws --endpoint-url $EP s3api head-object --bucket $B --key "$K" --query Metadata.md5chksum --output text | base64 -d | xxd -p)
check() { got=$(md5sum "$1" | cut -d' ' -f1); [[ "$want" == "$got" ]] && echo "   md5 OK" || echo "   md5 BAD"; rm -f "$1"; }
rate() { python3 -c "import os;b=os.path.getsize('$1');print(f'   {b/1e6:.0f} MB in {$3-$2:.0f}s = {b/1e6/($3-$2):.1f} MB/s')"; }
for n in $CONNS; do
  echo "s5cmd concurrency=$n part-size=8MB"
  t0=$(date +%s.%N)
  s5cmd --endpoint-url $EP --log error cp --concurrency $n --part-size 8 "s3://$B/$K" "$DEST/s5.gz"
  t1=$(date +%s.%N); rate "$DEST/s5.gz" $t0 $t1; check "$DEST/s5.gz"
  echo "rclone multi-thread-streams=$n chunk=8M"
  t0=$(date +%s.%N)
  rclone copyto ":s3,provider=Other,endpoint='$EP',env_auth=true:$B/$K" "$DEST/rc.gz" \
    --multi-thread-streams $n --multi-thread-chunk-size 8M --multi-thread-cutoff 16M --s3-no-check-bucket --ignore-checksum
  t1=$(date +%s.%N); rate "$DEST/rc.gz" $t0 $t1; check "$DEST/rc.gz"
done
