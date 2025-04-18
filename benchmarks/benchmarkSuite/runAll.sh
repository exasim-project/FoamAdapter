#!/bin/sh
cd "${0%/*}" || exit

current_dir=$(pwd)

cd implicitOperators || exit
snakemake -c1
cd "$current_dir" || exit

cd explicitOperators || exit
snakemake -c1
cd "$current_dir" || exit
