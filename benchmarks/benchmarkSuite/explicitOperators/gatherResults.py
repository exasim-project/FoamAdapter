# %%
from __future__ import annotations

from pathlib import Path
import pandas as pd

from foamlib.postprocessing.load_tables import datafile, load_tables
from foamlib.postprocessing.table_reader import read_catch2_benchmark

root = Path(__file__).parent
results = root / "results"
results.mkdir(exist_ok=True)

file = datafile(file_name="explicitOperators.xml", folder=".")
benchmark_results = load_tables(
    source=file, dir_name=root / "Cases", reader_fn=read_catch2_benchmark
)


# save per test case
def save_test_results(df: pd.DataFrame, test_case: str):
    test_case_df = df[df["test_case"] == test_case]
    if not test_case_df.empty:
        test_case_df.to_csv(results / f"{test_case}.csv", index=False)


for test_case in benchmark_results["test_case"].unique():
    save_test_results(benchmark_results, test_case)


# %%
