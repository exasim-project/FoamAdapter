from typing import Literal, Union, Annotated, Optional
from pydantic import Field, RootModel
from pybFoam.io.model_base import IOModelBase, IOModelMixin

# ---------- LES models (example) ----------
class Smagorinsky(IOModelBase):
    LESModel: Literal["Smagorinsky"]
    turbulence: Literal["on", "off"]
    delta: str = "cubeRootVol"
    Cs: float = 0.17


class OneEqEddy(IOModelBase):
    LESModel: Literal["oneEqEddy"]
    turbulence: Literal["on", "off"]
    delta: str = "cubeRootVol"
    Prt: float = 0.9

class LESConfig(
    RootModel[
        Annotated[
            Union[Smagorinsky, OneEqEddy],
            Field(discriminator="LESModel"),
        ]
    ],
    IOModelMixin
):
    @classmethod
    def from_ofdict(cls, d):
        """Parse LES config from OpenFOAM dictionary."""
        les_model_name = d.get[str]("LESModel")

        # Map model names to classes
        model_classes = {
            "Smagorinsky": Smagorinsky,
            "oneEqEddy": OneEqEddy,
        }
        
        if les_model_name not in model_classes:
            raise ValueError(f"Unknown LES model: {les_model_name}")
        
        model_class = model_classes[les_model_name]
        model_instance = model_class.from_ofdict(d)
        
        return cls(model_instance)