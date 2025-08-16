from dataclasses import dataclass
from typing import Annotated, Dict, Tuple, Type
from pydantic import BaseModel, Field, create_model

@dataclass(frozen=True)
class FileSpec:
    relpath: str
    encoding: str = "utf-8"
    required: bool = True

Registry = Dict[str, Tuple[Type[BaseModel], FileSpec]]


def build_case_inputs_class(
    reg: Registry,
    name: str = "CaseInputs",
) -> Type[BaseModel]:
    fields = {}
    for field_name, (model_type, spec) in reg.items():
        ann = Annotated[model_type, spec]                     # attach FileSpec via Annotated
        fields[field_name] = (ann, ...)                       # "..." means required
    CaseInputs = create_model(
        name,
        # __base__=base,
        **fields,
    )
    # CaseInputs.__module__ = __name__                          # improves pickling/docs
    return CaseInputs

