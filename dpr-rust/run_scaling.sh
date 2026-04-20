#!/bin/bash
# Finite-size scaling: extract S (order parameter jump, raw nodes) and C/N
# (largest cluster fraction at critical point) for logarithmically spaced N.
#
# N sweeps from 100 to 1e8 with one point per order of magnitude and one
# halfway between each (i.e. at multiples of sqrt(10) ~ 3.162).
#
# Rayon thread count is automatically capped at large N to stay within
# available RAM (~30 bytes * N * threads).

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BINARY="${1:-$SCRIPT_DIR/target/release/dpr}"
MEAN_CRIT=0.763
PROCESS=dpr
ALPHA=0.0
TOTAL_CPUS=$(nproc)

# Log-spaced N: every order of magnitude and halfway between (sqrt(10) steps)
N_VALUES=(
    100 316
    1000 3162
    10000 31623
    100000 316228
    1000000 3162278
    10000000 31622777
    100000000
)

# Run counts — generous since the user doesn't mind waiting.
runs_for_n() {
    local n=$1
    if   [ "$n" -ge 10000000 ]; then echo 10
    elif [ "$n" -ge  1000000 ]; then echo 20
    elif [ "$n" -ge   100000 ]; then echo 50
    elif [ "$n" -ge    10000 ]; then echo 100
    else                            echo 200
    fi
}

# Cap Rayon threads to keep peak RAM within ~75% of available memory.
# Each thread needs ~30 bytes * N; available RAM is read from /proc/meminfo.
threads_for_n() {
    local n=$1
    local ram_bytes
    ram_bytes=$(awk '/MemTotal/ {print $2 * 1024}' /proc/meminfo)
    local mem_per_thread=$(( n * 30 ))
    local max_threads=$(( ram_bytes * 3 / 4 / mem_per_thread ))
    # Clamp between 1 and total CPU count
    if   [ "$max_threads" -lt 1 ];              then echo 1
    elif [ "$max_threads" -gt "$TOTAL_CPUS" ];  then echo "$TOTAL_CPUS"
    else                                             echo "$max_threads"
    fi
}

OUTPUT=scaling.csv
echo "N,num_runs,threads,S,C_over_N,C_over_N_std,time_per_run_s" > "$OUTPUT"
echo "Detected $TOTAL_CPUS CPUs"
echo "Writing results to $OUTPUT"
echo "---"

for N in "${N_VALUES[@]}"; do
    NUM_RUNS=$(runs_for_n "$N")
    NUM_THREADS=$(threads_for_n "$N")

    # Window: 0.6*N to 0.9*N, bracketing the critical point at 0.763*N
    S_POINT=$(echo "$N * 6 / 10" | bc)
    E_POINT=$(echo "$N * 9 / 10" | bc)
    LEN=$E_POINT

    TMPFILE=$(mktemp /tmp/pm_${N}_XXXX.csv)

    echo -n "N=$N (runs=$NUM_RUNS, threads=$NUM_THREADS) ... "

    T_START=$(date +%s%N)

    RAYON_NUM_THREADS=$NUM_THREADS $BINARY perc-max \
        -N "$N" -k 2 \
        --s-point-f "$S_POINT" --e-point-f "$E_POINT" \
        --s-point-b "$S_POINT" --e-point-b "$E_POINT" \
        -r "$NUM_RUNS" -l "$LEN" \
        --mean-crit "$MEAN_CRIT" \
        -p "$PROCESS" -a "$ALPHA" \
        -o "$TMPFILE" 2>/dev/null

    T_END=$(date +%s%N)
    TIME_PER_RUN=$(awk "BEGIN {printf \"%.4f\", ($T_END - $T_START) / 1e9 / $NUM_RUNS}")

    NROWS=$(tail -n +2 "$TMPFILE" | wc -l)
    if [ "$NROWS" -lt 2 ]; then
        echo "ERROR: simulation produced $NROWS row(s), skipping"
        echo "$N,$NUM_RUNS,$NUM_THREADS,ERROR,ERROR,ERROR,ERROR" >> "$OUTPUT"
        rm -f "$TMPFILE"
        continue
    fi

    # jump_forward (col 2): ensemble mean S fraction — same on every row.
    # max_clust_mean_fld (col 6): per-run C/N — compute mean and std.
    STATS=$(awk -F',' -v N="$N" '
        NR == 1 { next }
        {
            s       = $2 * N
            cn      = $6 + 0
            sum_cn  += cn
            sum_cn2 += cn * cn
            n++
        }
        END {
            mean_cn = sum_cn / n
            std_cn  = (n > 1) ? sqrt((sum_cn2 - n * mean_cn^2) / (n - 1)) : 0
            print s, mean_cn, std_cn
        }
    ' "$TMPFILE")

    S_MEAN=$(echo $STATS  | awk '{print $1}')
    CN_MEAN=$(echo $STATS | awk '{print $2}')
    CN_STD=$(echo $STATS  | awk '{print $3}')

    echo "S(raw)=$S_MEAN  C/N=$CN_MEAN +/- $CN_STD  time/run=${TIME_PER_RUN}s"
    echo "$N,$NUM_RUNS,$NUM_THREADS,$S_MEAN,$CN_MEAN,$CN_STD,$TIME_PER_RUN" >> "$OUTPUT"

    rm -f "$TMPFILE"
done

echo "---"
echo "Done. Results in $OUTPUT"
