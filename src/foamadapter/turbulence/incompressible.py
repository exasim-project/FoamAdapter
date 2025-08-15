from typing import Literal, Union, Annotated, Optional
from pydantic import BaseModel, Field, RootModel
from pybFoam.io.model_base import IOModelBase
from foamadapter.turbulence.RASModel.incompressible import RASConfig
from foamadapter.turbulence.LESModel.incompressible import LESConfig



class LESProperties(IOModelBase):
    simulationType: Literal["LES"]
    LES: LESConfig


class RASProperties(IOModelBase):
    simulationType: Literal["RAS"]
    RAS: RASConfig


# ---------- Laminar ----------
class LaminarProperties(IOModelBase):
    simulationType: Literal["laminar"]
    # No additional block


# ---------- Single entry-point model with top-level discriminator ----------
class TurbulenceModel(
    RootModel[
        Annotated[
            Union[RASProperties, LESProperties, LaminarProperties],
            Field(discriminator="simulationType"),
        ]
    ]
):

    # Optional convenience accessors:
    @property
    def simulationType(self) -> str:
        return self.__root__.simulationType

    @property
    def RAS(self) -> Optional[RASConfig]:
        return getattr(self.__root__, "RAS", None)

    @property
    def LES(self) -> Optional[LESConfig]:
        return getattr(self.__root__, "LES", None)


from pprint import pprint

# Generate JSON Schema
schema = TurbulenceModel.model_json_schema()

# Pretty-print
pprint(schema)

import json

with open("turbulence_schema.json", "w") as f:
    json.dump(schema, f, indent=2)
