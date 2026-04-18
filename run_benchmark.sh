#!/bin/bash
# Benchmark the DPR binary across a range of N values.
# Reports wall-clock time per single run (all CPUs used), edge throughput
# (edges/second), and memory per thread estimate.
#
# Usage:
#   ./run_benchmark.sh [binary]
#
# Example:
#   ./run_benchmark.sh ./target/release/dpr

BINARY="${1:-./target/release/dpr}"
TOTAL_CPUS=$(nproc)
MEAN_CRIT=0.763

if [ ! -x "$BINARY" ]; then
    echo "ERROR: binary not found or not executable: $BINARY"
    echo "Build with: cargo build --release"
    exit 1
fi

# N values: one per order of magnitude for a quick sweep.
N_VALUES=(1000 10000 100000 1000000 10000000)

# Fixed run count per N — enough for a stable mean without taking too long.
RUNS_PER_N=10

# Cap threads the same way run_scaling.sh does (~30 bytes * N per thread).
threads_for_n() {
    local n=$1
    local ram_bytes
    ram_bytes=$(awk '/MemTotal/ {print $2 * 1024}' /proc/meminfo)
    local mem_per_thread=$(( n * 30 ))
    local max_threads=$(( ram_bytes * 3 / 4 / mem_per_thread ))
    if   [ "$max_threads" -lt 1 ];             then echo 1
    elif [ "$max_threads" -gt "$TOTAL_CPUS" ]; then echo "$TOTAL_CPUS"
    else                                            echo "$max_threads"
    fi
}

OUTPUT=benchmark.csv
echo "N,threads,runs,time_per_run_s,edges_per_s,mem_per_thread_MB" > "$OUTPUT"

echo "DPR benchmark  —  binary: $BINARY"
echo "CPUs available: $TOTAL_CPUS"
printf "%-12s  %-8s  %-14s  %-18s  %-18s\n" \
    "N" "threads" "time/run (s)" "edges/s (M)" "mem/thread (MB)"
echo "------------------------------------------------------------------------"

for N in "${N_VALUES[@]}"; do
    NUM_THREADS=$(threads_for_n "$N")

    # Simulate 0.9 * N edges to cover the critical point window.
    LEN=$(echo "$N * 9 / 10" | bc)
    S_POINT=$(echo "$N * 6 / 10" | bc)
    E_POINT=$LEN

    TMPFILE=$(mktemp /tmp/bench_${N}_XXXX.csv)

    T_START=$(date +%s%N)

    RAYON_NUM_THREADS=$NUM_THREADS "$BINARY" perc-max \
        -N "$N" -k 2 \
        --s-point-f "$S_POINT" --e-point-f "$E_POINT" \
        --s-point-b "$S_POINT" --e-point-b "$E_POINT" \
        -r "$RUNS_PER_N" -l "$LEN" \
        --mean-crit "$MEAN_CRIT" \
        -p dpr -a 0.0 \
        -o "$TMPFILE" >/dev/null 2>&1

    EXIT_CODE=$?
    T_END=$(date +%s%N)

    rm -f "$TMPFILE"

    if [ "$EXIT_CODE" -ne 0 ]; then
        echo "N=$N  ERROR (exit $EXIT_CODE)"
        echo "$N,$NUM_THREADS,$RUNS_PER_N,ERROR,ERROR,ERROR" >> "$OUTPUT"
        continue
    fi

    # Wall time is for all runs in parallel; divide by run count for per-run time.
    TIME_TOTAL_S=$(awk "BEGIN {printf \"%.6f\", ($T_END - $T_START) / 1e9}")
    TIME_PER_RUN=$(awk "BEGIN {printf \"%.4f\", $TIME_TOTAL_S / $RUNS_PER_N}")

    # Throughput: edges simulated per second across the whole batch.
    # Each run processes LEN forward + LEN backward edges = 2 * LEN.
    TOTAL_EDGES=$(( RUNS_PER_N * LEN * 2 ))
    EDGES_PER_S=$(awk "BEGIN {printf \"%.2f\", $TOTAL_EDGES / $TIME_TOTAL_S / 1e6}")

    # Rough memory per thread: degree array (4B * N) + edge arrays (8B * LEN)
    # + union-find (9B * N) + edge set if N < 50000 (14B * LEN).
    MEM_MB=$(awk -v n="$N" -v l="$LEN" 'BEGIN {
        uf   = 9 * n
        deg  = 4 * n
        edge = 8 * l
        dedup = (n < 50000) ? 14 * l : 0
        printf "%.1f", (uf + deg + edge + dedup) / 1048576
    }')

    printf "%-12s  %-8s  %-14s  %-18s  %-18s\n" \
        "$N" "$NUM_THREADS" "${TIME_PER_RUN}s" "${EDGES_PER_S}M/s" "${MEM_MB} MB"

    echo "$N,$NUM_THREADS,$RUNS_PER_N,$TIME_PER_RUN,$EDGES_PER_S,$MEM_MB" >> "$OUTPUT"
done

echo "------------------------------------------------------------------------"
echo "Results saved to $OUTPUT"
