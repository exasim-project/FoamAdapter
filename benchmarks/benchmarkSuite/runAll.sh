#!/bin/bash
cd "${0%/*}" || exit

current_dir=$(pwd)

# run benchmarks
run_benchmark() {
    echo "Running $2 benchmarks in: $1"
    # mkdir -p $1/results
    find $1 -name "Allrun" -exec {} bench_$2 \;

    echo "Gathering results..." $1
    python3 gatherResults.py $1

    echo "Completed benchmark: $1"
}

# Execute the benchmark commands
echo "Creating study..."
python3 createStudies.py

# Define benchmarks to run
benchmarks=("explicitOperators" "implicitOperators" "dsl")
# Run each benchmark
for benchmark in "${benchmarks[@]}"; do
    run_benchmark "$current_dir/$benchmark" $benchmark
done
