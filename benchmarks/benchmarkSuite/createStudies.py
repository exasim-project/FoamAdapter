# %%
# create 2DSquare and 3DCube benchmark cases

from pathlib import Path
from foamlib.preprocessing.parameter_study import record_generator


def build_records(name, resolution):
    return [ {"case_name": f"{name}_N{str(r)}", "Res":r, "MeshType": name, "Resolution": f"N{r}" } for r in resolution ]

def create_cases(root):
    study_cube = record_generator(
        records=build_records("3DCube", [10, 20, 50, 100, 200]),
        template_case=root / "templates/3DCube",
        output_folder=root / "Cases",
    )

    template_case = root / "templates/2DSquare"
    study_square = record_generator(
        records=build_records("2DSquare", [10, 20, 50, 100, 200, 500, 1000, 2000]),
        template_case=template_case,
        output_folder=root / "Cases"
    )

    # Combine both studies
    final_study = study_cube + study_square
    final_study.create_study(study_base_folder=root)

root = Path(__file__).parent
for c in ["dsl",  "explicitOperators",  "implicitOperators"]:
    create_cases(root/c)
