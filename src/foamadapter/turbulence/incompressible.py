from typing import Literal, Union, Annotated, Optional
from pydantic import BaseModel, Field, RootModel
from pybFoam.io.model_base import IOModelMixin, IOModelBase
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
    ],
    IOModelMixin
):
    
    @classmethod
    def from_ofdict(cls, d):
        """Override from_ofdict for RootModel compatibility."""
        # For RootModel, we need to parse the entire dictionary content
        # instead of trying to extract named fields
        
        # First, get the simulationType to determine which model to use
        simulation_type = d.get[str]("simulationType")
        
        # Build the mapping for the specific model type
        mapping = {"simulationType": simulation_type}
        
        # Add type-specific fields
        if simulation_type == "RAS" and d.isDict("RAS"):
            from foamadapter.turbulence.RASModel.incompressible import RASConfig
            mapping["RAS"] = RASConfig.from_ofdict(d.subDict("RAS"))
        elif simulation_type == "LES" and d.isDict("LES"):
            from foamadapter.turbulence.LESModel.incompressible import LESConfig
            mapping["LES"] = LESConfig.from_ofdict(d.subDict("LES"))
        # For laminar, no additional fields needed
        
        return cls(mapping)

    # # Optional convenience accessors:
    @property
    def simulationType(self) -> str:
        return self.root.simulationType

    @property
    def RAS(self) -> Optional[RASConfig]:
        return getattr(self.root, "RAS", None)

    @property
    def LES(self) -> Optional[LESConfig]:
        return getattr(self.root, "LES", None)


# from pprint import pprint

# # Generate JSON Schema
# schema = TurbulenceModel.model_json_schema()

# # Pretty-print
# pprint(schema)

# import json

# with open("turbulence_schema.json", "w") as f:
#     json.dump(schema, f, indent=2)
