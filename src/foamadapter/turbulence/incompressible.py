from typing import Literal, Union, Annotated, Optional
from pydantic import BaseModel, Field, RootModel
from pybFoam.io.model_base import IOModelMixin, IOModelBase
from pybFoam import dictionary
from pybFoam.turbulence import incompressibleTurbulenceModel
from foamadapter.turbulence.RASModel.incompressible import RASConfig
from foamadapter.turbulence.LESModel.incompressible import LESConfig
from foamadapter.inputs_files.case_inputs import Registry, FileSpec

from foamadapter.modules import fields
from foamadapter.modules import models



class LESProperties(IOModelBase):
    simulationType: Literal["LES"]
    LES: LESConfig

    @classmethod
    def additional_inputs(cls, reg: Registry) -> Registry:
        """Example method for LES properties."""
        reg["LES"] = 0
        return reg


class RASProperties(IOModelBase):
    simulationType: Literal["RAS"]
    RAS: RASConfig

    @classmethod
    def additional_inputs(cls, reg: Registry) -> Registry:
        """Example method for RAS properties."""
        reg["RAS"] = 1
        return reg


class LaminarProperties(IOModelBase):
    simulationType: Literal["laminar"]
    # No additional block

    @classmethod
    def additional_inputs(cls, reg: Registry) -> Registry:
        """Example method for laminar properties."""
        return reg


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
    
    @property
    def dependencies(self) -> list[str]:
        return ["U", "phi", "singlePhaseTransportModel"]
    
    @property
    def description(self) -> str:
        return "Incompressible turbulence model"
    
    def __call__(self, deps: dict):
        U = fields.get_field(deps, "U")
        phi = fields.get_field(deps, "phi")
        singlePhaseTransportModel = models.get_model(deps, "singlePhaseTransportModel")
        turbulence = incompressibleTurbulenceModel.New(U, phi, singlePhaseTransportModel)
        return turbulence


    @classmethod
    def inputs(cls, reg: Registry, turb_file: str = "constant/turbulenceProperties") -> Registry:
        """Return the input parameters for the turbulence model."""
        reg["turbulenceProperties"] = (
            cls,
            FileSpec(turb_file, encoding="utf-8", required=True)
        )

        d = dictionary.read(turb_file)

        simulation_type = d.get[str]("simulationType")
      
        # Add type-specific fields
        if simulation_type == "RAS" and d.isDict("RAS"):
            from foamadapter.turbulence.RASModel.incompressible import RASConfig
            RASConfig.additional_inputs(d.subDict("RAS"), reg)
        elif simulation_type == "LES" and d.isDict("LES"):
            from foamadapter.turbulence.LESModel.incompressible import LESConfig
            LESConfig.additional_inputs(d.subDict("LES"), reg)

        return reg