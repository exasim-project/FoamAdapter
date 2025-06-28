# %%

from pathlib import Path

from foamlib.preprocessing.parameter_study import csv_generator


root = Path(__file__).parent

# 3DCube
template_case = root / "templates/3DCube"
study_cube = csv_generator(
    csv_file=root / "parastudy_3DCube.csv", template_case=template_case, output_folder=root / "Cases"
)

# 2DSquare
template_case = root / "templates/2DSquare"
study_square = csv_generator(
    csv_file=root / "parastudy_2DSquare.csv", template_case=template_case, output_folder=root / "Cases"
)

# Combine both studies
final_study = study_cube + study_square
final_study.create_study(study_base_folder=root)

