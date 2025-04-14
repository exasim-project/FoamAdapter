# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the CSV file
df = pd.read_csv("results/benchmark_results.csv")

# Create a label for plotting
df["Label"] = df["Geometry"] + "_" + df["Resolution"] + " | " + df["Executor"]


def normalize_all_by_openfoam(group):
    openfoam_rows = group[group["Executor"] == "OpenFOAM"]
    if not openfoam_rows.empty:
        openfoam_mean = openfoam_rows["Mean"].iloc[0]
        group["Normalized_Mean"] = group["Mean"] / openfoam_mean
    else:
        group["Normalized_Mean"] = float("nan")
    return group


# Apply corrected normalization across Geometry + Resolution
df_fixed_normalization = df.groupby(["TestCase","Geometry", "Resolution"], group_keys=False).apply(
    normalize_all_by_openfoam
)
#%%
# Plot using seaborn
df_fixed_normalization["Label"] = (
    df_fixed_normalization["Geometry"]
    + "_"
    + df_fixed_normalization["Resolution"]
    + " | "
    + df_fixed_normalization["Executor"]
)
df_fixed_normalization["Speed up"] = 1 / df_fixed_normalization["Normalized_Mean"]
# resolution_order = ["N10", "N20", "N50", "N100", "N200", "N500", "N1000", "N1000"]
# df_fixed_normalization["Resolution"] = pd.Categorical(
#     df_fixed_normalization["Resolution"],
#     categories=resolution_order,
#     ordered=True
# )
# %%
# Plot 3D cube
def plot(df_fixed_normalization, testcase):
    # Filter the DataFrame for the specific test case
    df = df_fixed_normalization[df_fixed_normalization["TestCase"] == testcase]
    df_fixed_normalization_3dCube = df[df["Geometry"] == "3DCube"]
    resolution_order = ["N10", "N20", "N50", "N100", "N200"]
    df_fixed_normalization_3dCube["Resolution"] = pd.Categorical(
        df_fixed_normalization_3dCube["Resolution"],
        categories=resolution_order,
        ordered=True
    )

    plt.figure(figsize=(14, 6))
    sns.barplot(
        data=df_fixed_normalization_3dCube,
        x="Resolution",
        y="Speed up",
        hue="Executor"
    )
    # plt.xticks(rotation=90)
    plt.axhline(1, color='black', linestyle='--', linewidth=1)  # Add baseline
    plt.ylabel("Normalized Mean (w.r.t. OpenFOAM = 1)")
    plt.title(f"3DCube: Speed up for the explicit {testcase}")
    plt.tight_layout()
    plt.savefig("results/{testcase}_3DCube.png")

    df_fixed_normalization_2DSquare = df[df["Geometry"] == "2DSquare"]
    resolution_order = ["N20", "N50", "N100", "N200", "N500", "N1000", "N2000"]
    df_fixed_normalization_2DSquare["Resolution"] = pd.Categorical(
        df_fixed_normalization_2DSquare["Resolution"],
        categories=resolution_order,
        ordered=True
    )
    plt.figure(figsize=(14, 6))
    sns.barplot(
        data=df_fixed_normalization_2DSquare,
        x="Resolution",
        y="Speed up",
        hue="Executor"
    )
    # plt.xticks(rotation=90)
    plt.axhline(1, color='black', linestyle='--', linewidth=1)  # Add baseline
    plt.ylabel("Normalized Mean (w.r.t. OpenFOAM = 1)")
    plt.title(f"2DSquare: Speed up for the explicit {testcase}")
    plt.tight_layout()
    plt.savefig(f"results/{testcase}_2DSquare.png")
    plt.show()
# %%
plot(df_fixed_normalization, "GradOperator")
plot(df_fixed_normalization, "DivOperator")
# %%
