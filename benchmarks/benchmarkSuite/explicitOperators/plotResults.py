# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# %%

## plot results as bar plots


def create_bar_plot(table, testcase, mesh_type, resolution_order):
    sliced_table = table[table["MeshType"] == mesh_type].copy()
    sliced_table["Resolution"] = pd.Categorical(
        sliced_table.loc[:, "Resolution"],
        categories=resolution_order,
        ordered=True
    )
    plt.figure(figsize=(14, 6))
    sliced_table.loc[:, "Label"] = (
        sliced_table.loc[:, "benchmark_name"] + "_" + sliced_table.loc[:, "section2"]
    )
    sns.barplot(
        data=sliced_table,
        x="Resolution",
        y="normalized_speedup",
        hue="Label"
    )
    # plt.xticks(rotation=90)
    plt.axhline(1, color="black", linestyle="--", linewidth=1)  # Add baseline
    plt.ylabel("Normalized Mean (w.r.t. OpenFOAM = 1)")
    plt.title(f"{mesh_type}: Speed up for the explicit {testcase}")
    plt.tight_layout()
    plt.yscale("log")
    plt.savefig(f"results/{testcase}_{mesh_type}.png")


def plot_table(test_case):
    # Load the CSV file
    table = pd.read_csv(Path(f"results/{test_case}.csv"))

    create_bar_plot(
        table=table,
        testcase=test_case,
        mesh_type="3DCube",
        resolution_order=["N10", "N20", "N50", "N100", "N200"],
    )

    create_bar_plot(
        table=table,
        testcase=test_case,
        mesh_type="2DSquare",
        resolution_order=["N20", "N50", "N100", "N200", "N500", "N1000", "N2000"],
    )


plot_table("DivOperator")
plot_table("LaplacianOperator")
plot_table("GradOperator")
plot_table("FaceInterpolation")
plot_table("FaceNormalGradient")
