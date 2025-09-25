# %%
from __future__ import annotations

from pathlib import Path
import os

from foamlib.postprocessing.load_tables import datafile, load_tables
from foamlib.postprocessing.table_reader import read_catch2_benchmark

# Function to normalize avg_runtime within each group relative to OpenFOAM
def normalize_group(group):
    baseline = group.loc[group["benchmark_name"] == "OpenFOAM", "avg_runtime"]
    if not baseline.empty:
        group["normalized_speedup"] = baseline.values[0] / group["avg_runtime"]
    else:
        group["normalized_speedup"] = float("nan")  # No baseline found
    return group


# save per test case
def save_test_results(df, test_case: str):
    test_case_df = df[df["test_case"] == test_case]
    if not test_case_df.empty:
        test_case_df = test_case_df.groupby(group_keys).apply(
            normalize_group, include_groups=False
        ).reset_index()
        test_case_df.to_csv(results / f"{test_case}.csv", index=False)

if __name__ == "__main__":
    root = Path(os.getcwd())
    results = root / "results"
    results.mkdir(exist_ok=True)

    cases = root / "Cases"
    if not cases.exists():
        print(f"could not find {cases}")

    group_keys = ["MeshType", "Resolution"]

    file = datafile(file_name="result.xml", folder=".")
    benchmark_results = load_tables(
        source=file, dir_name=cases, reader_fn=read_catch2_benchmark
    )
    if benchmark_results:
        for test_case in benchmark_results["test_case"].unique():
            save_test_results(benchmark_results, test_case)
    else:
        print(f"Failed to postprocess {cases}")



# %%
