#!/bin/bash
cd "${0%/*}" || exit

current_dir=$(pwd)

# remove Cases directories
# Define benchmarks to run
benchmarks=("dsl" "explicitOperators" "implicitOperators")

# Run each benchmark
for benchmark in "${benchmarks[@]}"; do
    rm -rf "$current_dir/$benchmark"
done
