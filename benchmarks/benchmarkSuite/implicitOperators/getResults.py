#%%
# Re-run code after environment reset
import pandas as pd
import os
from pathlib import Path
import xml.etree.ElementTree as ET

def parse_catch2_benchmark_xml(xml_path, geometry_label, resolution_label):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    rows = []

    for testcase in root.findall('TestCase'):
        test_name = testcase.get('name')
        for section in testcase.findall('Section'):
            section_name = section.get('name')
            for sub_section in section.findall('Section'):
                sub_section_name = sub_section.get('name')
                results = sub_section.find('BenchmarkResults')
                if results is None:
                    continue

                mean = results.find('mean')
                stddev = results.find('standardDeviation')
                outliers = results.find('outliers')

                row = {
                    'TestCase': test_name,
                    'Section': section_name,
                    'SubSection': sub_section_name,
                    'Executor': results.get('name'),
                    'Mean': float(mean.get('value')),
                    # 'Mean_Lower': float(mean.get('lowerBound')),
                    # 'Mean_Upper': float(mean.get('upperBound')),
                    # 'StdDev': float(stddev.get('value')),
                    # 'StdDev_Lower': float(stddev.get('lowerBound')),
                    # 'StdDev_Upper': float(stddev.get('upperBound')),
                    # 'Outliers_HighMild': int(outliers.get('highMild')),
                    # 'Outliers_HighSevere': int(outliers.get('highSevere')),
                    'Geometry': geometry_label,
                    'Resolution': resolution_label
                }
                rows.append(row)
    return pd.DataFrame(rows)

# %%
all_dfs = []
basedir = "cases/implicitOperators"
for folder in os.listdir(basedir):
    print("Processing folder:", folder)
    xml_path = f"{basedir}/{folder}/implicitOperators.xml"
    geometry = folder.split("_")[0]
    resolution = folder.split("_")[-1]  # e.g., N10, N20
    df = parse_catch2_benchmark_xml(xml_path, geometry, resolution)
    all_dfs.append(df)

combined_df = pd.concat(all_dfs, ignore_index=True)
# %%
Path("results").mkdir(exist_ok=True)
combined_df.to_csv("results/benchmark_results.csv", index=False)
