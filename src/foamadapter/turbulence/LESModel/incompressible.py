from typing import Literal, Union, Annotated, Optional
from pydantic import Field
from pybFoam.io.model_base import IOModelBase

# ---------- LES models (example) ----------
class LES_Smagorinsky(IOModelBase):
    LESModel: Literal["Smagorinsky"]
    turbulence: Literal["on", "off"]
    delta: str = "cubeRootVol"
    Cs: float = 0.17


class LES_OneEqEddy(IOModelBase):
    LESModel: Literal["oneEqEddy"]
    turbulence: Literal["on", "off"]
    delta: str = "cubeRootVol"
    Prt: float = 0.9


LESConfig = Annotated[
    Union[LES_Smagorinsky, LES_OneEqEddy],
    Field(discriminator="LESModel"),
]


