#!/bin/bash
# Finite-size scaling: extract S (jump size) and C/N (largest cluster at critical
# point) for varying N, then write results to scaling.csv.
# S is expressed as raw node count (fraction * N).

BINARY="${1:-./target/release/dpr}"
MEAN_CRIT=0.763
PROCESS=dpr
ALPHA=0.0

# N values to sweep — add or remove as needed
N_VALUES=(500 1000 2000 5000 10000 20000 50000 100000 500000 1000000 5000000)

# Scale runs down for large N to keep runtime manageable.
# Adjust the thresholds to taste.
runs_for_n() {
    local n=$1
    if   [ "$n" -ge 1000000 ]; then echo 20
    elif [ "$n" -ge  100000 ]; then echo 50
    elif [ "$n" -ge   10000 ]; then echo 100
    else                            echo 200
    fi
}

OUTPUT=scaling.csv
echo "N,num_runs,S,S_std,C_over_N,C_over_N_std" > "$OUTPUT"
echo "Writing results to $OUTPUT"
echo "---"

for N in "${N_VALUES[@]}"; do
    NUM_RUNS=$(runs_for_n "$N")

    # Window spans roughly 0.6*N to 0.9*N, centred on critical point 0.763*N
    S_POINT=$(echo "$N * 6 / 10" | bc)
    E_POINT=$(echo "$N * 9 / 10" | bc)
    LEN=$E_POINT

    TMPFILE=$(mktemp /tmp/pm_${N}_XXXX.csv)

    echo -n "N=$N (runs=$NUM_RUNS) ... "

    $BINARY perc-max \
        -N "$N" -k 2 \
        --s-point-f "$S_POINT" --e-point-f "$E_POINT" \
        --s-point-b "$S_POINT" --e-point-b "$E_POINT" \
        -r "$NUM_RUNS" -l "$LEN" \
        --mean-crit "$MEAN_CRIT" \
        -p "$PROCESS" -a "$ALPHA" \
        -o "$TMPFILE" 2>/dev/null

    # Check the simulation produced output before parsing.
    NROWS=$(tail -n +2 "$TMPFILE" | wc -l)
    if [ "$NROWS" -lt 2 ]; then
        echo "ERROR: simulation produced $NROWS row(s), skipping"
        echo "$N,$NUM_RUNS,ERROR,ERROR,ERROR,ERROR" >> "$OUTPUT"
        rm -f "$TMPFILE"
        continue
    fi

    # jump_forward (col 2) is the ensemble mean S — same on every data row.
    # max_clust_mean_fld (col 6) is per-run — compute mean and std across runs.
    STATS=$(awk -F',' -v N="$N" '
        NR == 1 { next }
        {
            s       = $2 * N    # S as raw node count
            cn      = $6 + 0
            sum_cn  += cn
            sum_cn2 += cn * cn
            n++
        }
        END {
            if (n < 2) { print s, 0, sum_cn/n, 0; exit }
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
    echo "$N,$NUM_RUNS,$S_MEAN,$S_STD,$CN_MEAN,$CN_STD" >> "$OUTPUT"

    rm -f "$TMPFILE"
done

echo "---"
echo "Done. Results in $OUTPUT"
