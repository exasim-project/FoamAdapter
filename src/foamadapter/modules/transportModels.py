from typing import Literal
from pydantic import Field
from pybFoam.io.model_base import IOModelBase
from pybFoam.turbulence import singlePhaseTransportModel
from foamadapter.inputs_files.case_inputs import Registry, FileSpec
from foamadapter.modules.models import Models, Model
from foamadapter.modules import fields


@Models.deps("U", "phi")
class SinglePhaseTransportModel(
    IOModelBase
):
    # nu: float = Field(..., description="Kinematic viscosity", gt=0)
    # transportModel: Literal["Newtonian"]

    @classmethod
    def from_ofdict(cls, d):
        """Override from_ofdict for RootModel compatibility."""
        # For RootModel, we need to parse the entire dictionary content
        # instead of trying to extract named fields
        return super().from_ofdict(d)

    
    def __call__(self, deps: dict):
        U = fields.get_field(deps, "U")
        phi = fields.get_field(deps, "phi")
        transportModel = singlePhaseTransportModel(U, phi)
        return transportModel

    @classmethod
    def inputs(cls, reg: Registry, transport_properties: str = "constant/transportProperties") -> Registry:
        """Return the input parameters for the turbulence model."""
        reg["transportProperties"] = (
            cls,
            FileSpec(transport_properties, encoding="utf-8", required=True)
        )

        return reg