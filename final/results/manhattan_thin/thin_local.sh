#!/bin/bash
# Thin one raw GWAS scan for Manhattan plotting.
#   usage: thin_local.sh <trait> <MODEL> <in_dir> <out_dir>
# Columns: 1 SNP 2 Chr 3 Pos 4 P.value ... 10 model
# Keeps every SNP with P < 0.01 plus a 1-in-25 subsample of the rest.
set -euo pipefail
TRAIT="$1"; MODEL="$2"; IN="$3"; OUT="$4"
mkdir -p "$OUT"
awk -F, -v m="$MODEL" '
  NR==1 { next }
  $10 != m { next }
  { n++
    p = $4 + 0
    if (p <= 0) next
    if ($2+0 < 1 || $2+0 > 10) next
    if (p < 0.01 || n % 25 == 0) { print $2","$3","p; kept++ }
  }
  END { printf "  %s/%s: %d rows -> %d kept\n", "'"$TRAIT"'", m, n, kept > "/dev/stderr" }
' "$IN/$TRAIT.csv" > "$OUT/${TRAIT}__${MODEL}.csv.tmp"
{ echo "Chr,Pos,P"; cat "$OUT/${TRAIT}__${MODEL}.csv.tmp"; } > "$OUT/${TRAIT}__${MODEL}.csv"
mv "$OUT/${TRAIT}__${MODEL}.csv.tmp" "$OUT/.tmp_${TRAIT}_${MODEL}" 2>/dev/null || true
