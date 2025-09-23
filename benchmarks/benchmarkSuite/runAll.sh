#!/bin/bash
cd "${0%/*}" || exit

current_dir=$(pwd)

# run benchmarks
run_benchmark() {
    local benchmark_dir="$1"
    echo "Running benchmark: $benchmark_dir"

    # Save the current directory
    local start_dir="$current_dir"

    # Navigate to the benchmark directory
    cd "$benchmark_dir" || { echo "Failed to cd to $benchmark_dir"; return 1; }

    echo "Running benchmark..."
    ./runAll.sh

    echo "Gathering results..." $1
    python3 ../gatherResults.py

    # Return to the original directory
    cd "$start_dir" || { echo "Failed to return to $start_dir"; return 1; }

    echo "Completed benchmark: $benchmark_dir"
}

# Execute the benchmark commands
echo "Creating study..."
python3 createStudies.py

# Define benchmarks to run
benchmarks=("explicitOperators" "implicitOperators" "dsl")
# Run each benchmark
for benchmark in "${benchmarks[@]}"; do
    run_benchmark "$current_dir/$benchmark"
done
