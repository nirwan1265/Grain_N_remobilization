#!/bin/bash
################################################################################
## Scan every split_by_trait CSV and report, per trait and per model:
##   total rows, hits at the suggestive threshold, hits at Bonferroni,
##   and the best (smallest) p-value that model achieved.
##
## The "best P" column is the one that answers "did this model have nothing, or
## was it never run?" -- a model that ran but found nothing still has a best P
## (somewhere around 1e-5); a model that never ran has no rows at all.
##
## Column layout of the input files (comma separated):
##   1 SNP  2 Chr  3 Pos  4 P.value  5 MAF  6 nobs  7 H.B.P.Value  8 Effect
##   9 trait  10 model  11 FDR  12 p_star  13 line_log10
##
## Usage:
##   bash check_models_and_hits.sh                      # all traits
##   bash check_models_and_hits.sh Q.E.csv R.csv        # just these
##
## One full pass per file over ~8-12M rows; budget ~30-60 s per trait, so the
## whole directory is roughly 40-60 min. Run it under bsub or nohup.
################################################################################

set -uo pipefail

IN_DIR="/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/split_by_trait"
M=4177796                       # markers tested
SUGG=$(awk -v m=$M 'BEGIN{printf "%.10e", 1/m}')
BONF=$(awk -v m=$M 'BEGIN{printf "%.10e", 0.05/m}')

echo "# m=$M  suggestive(1/m)=$SUGG  bonferroni(0.05/m)=$BONF"
printf "%-14s %-8s %10s %8s %8s %12s %8s\n" \
       TRAIT MODEL ROWS SUGG BONF BEST_P NEGLOG10

cd "$IN_DIR" || exit 1
FILES=( "${@:-}" )
if [[ -z "${FILES[0]:-}" ]]; then FILES=( *.csv ); fi

for f in "${FILES[@]}"; do
  [[ -f "$f" ]] || { echo "MISSING FILE: $f" >&2; continue; }
  trait="${f%.csv}"
  awk -F, -v t="$trait" -v s="$SUGG" -v b="$BONF" '
    NR == 1 { next }
    {
      m = $10; p = $4 + 0
      n[m]++
      if (p > 0 && (!(m in best) || p < best[m])) best[m] = p
      if (p > 0 && p <= s) { sg[m]++; if (p <= b) bf[m]++ }
    }
    END {
      for (k in n)
        printf "%-14s %-8s %10d %8d %8d %12.3e %8.2f\n",
               t, k, n[k], sg[k] + 0, bf[k] + 0, best[k],
               -log(best[k]) / log(10)
    }' "$f" | sort -k2,2
done
