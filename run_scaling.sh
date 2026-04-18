#!/bin/bash
# Finite-size scaling: extract S (jump size) and C/N (largest cluster at critical
# point) for varying N, then write results to scaling.csv.

BINARY="${1:-./target/release/dpr}"
MEAN_CRIT=0.763
NUM_RUNS=200
PROCESS=dpr
ALPHA=0.0

# N values to sweep — add or remove as needed
N_VALUES=(500 1000 2000 5000 10000 20000)

OUTPUT=scaling.csv
echo "N,S,S_std,C_over_N,C_over_N_std" > "$OUTPUT"
# Note: S is expressed as raw node count (fraction * N)
echo "Writing results to $OUTPUT"
echo "---"

for N in "${N_VALUES[@]}"; do
    # Window spans roughly 0.6*N to 0.9*N, centred on critical point 0.763*N
    S_POINT=$(echo "$N * 6 / 10" | bc)
    E_POINT=$(echo "$N * 9 / 10" | bc)
    LEN=$E_POINT

    TMPFILE=$(mktemp /tmp/pm_${N}_XXXX.csv)

    echo -n "N=$N ... "

    $BINARY perc-max \
        -N "$N" -k 2 \
        --s-point-f "$S_POINT" --e-point-f "$E_POINT" \
        --s-point-b "$S_POINT" --e-point-b "$E_POINT" \
        -r "$NUM_RUNS" -l "$LEN" \
        --mean-crit "$MEAN_CRIT" \
        -p "$PROCESS" -a "$ALPHA" \
        -o "$TMPFILE" 2>/dev/null

    # jump_forward (col 2) is the ensemble mean S — same on every data row
    # jump_std would need re-adding to CSV; compute from per-run data instead.
    # max_clust_mean_fld (col 6) is per-run — average + std across runs.
    STATS=$(awk -F',' -v N="$N" '
        NR == 1 { next }   # skip header
        {
            s    = $2 * N   # multiply fraction by N for raw node count
            cn   = $6 + 0
            sum_cn  += cn
            sum_cn2 += cn * cn
            n++
        }
        END {
            mean_cn = sum_cn / n
            std_cn  = sqrt((sum_cn2 - n * mean_cn^2) / (n - 1))
            print s, 0, mean_cn, std_cn
        }
    ' "$TMPFILE")

    S_MEAN=$(echo $STATS | awk '{print $1}')
    S_STD=$(echo $STATS  | awk '{print $2}')
    CN_MEAN=$(echo $STATS | awk '{print $3}')
    CN_STD=$(echo $STATS  | awk '{print $4}')

    echo "S(raw)=$S_MEAN  C/N=$CN_MEAN"
    echo "$N,$S_MEAN,$S_STD,$CN_MEAN,$CN_STD" >> "$OUTPUT"

    rm -f "$TMPFILE"
done

echo "---"
echo "Done. Results in $OUTPUT"
