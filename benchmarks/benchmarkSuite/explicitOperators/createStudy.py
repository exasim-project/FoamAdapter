# %%

from pathlib import Path

from foamlib.preprocessing.parameter_study import csv_generator

# baseCase
root = Path(__file__).parent
template_case = root / "baseCase"

study = csv_generator(
    csv_file=root / "parastudy.csv", template_case=template_case, output_folder=root / "Cases"
)

study.create_study(study_base_folder=root)
