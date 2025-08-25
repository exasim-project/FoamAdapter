#!/bin/bash -l
set -eo pipefail

echo "Sourcing OpenFOAM..."
source /usr/local/app/openfoam/OpenFOAM-v2412/etc/bashrc

# -------------------------
# Step 1: Prepare NeoN
# -------------------------
echo "Cloning NeoN (main branch)..."
git clone --depth 1 --branch main https://gitlab-ce.lrz.de/greole/neon.git ../NeoN

# -------------------------
# Step 2: Build FoamAdapter
# -------------------------
echo "Building FoamAdapter against NeoN..."
cmake --preset develop -DFOAMADAPTER_NEON_DIR=../NeoN \
      -DCMAKE_CUDA_ARCHITECTURES=89 -DNeoN_WITH_THREADS=OFF
cmake --build --preset develop

# -------------------------
# Step 3: Run Tests
# -------------------------
echo "Running FoamAdapter tests..."
ctest --preset develop --output-on-failure
