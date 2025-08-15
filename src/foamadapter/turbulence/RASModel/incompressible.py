from typing import Literal, Union, Annotated, Optional
from pydantic import Field, RootModel
from pybFoam.io.model_base import IOModelBase, IOModelMixin


# ---------- RAS models (each with its own extra fields) ----------
class kEpsilon(IOModelBase):
    RASModel: Literal["kEpsilon"] = "kEpsilon"
    turbulence: Literal["on", "off"] = "on"
    printCoeffs: Literal["on", "off"] = "on"
    Cmu: float = 0.09
    C1: float = 1.44
    C2: float = 1.92
    sigma_k: float = 1.0
    sigma_epsilon: float = 1.3


class kOmega(IOModelBase):
    RASModel: Literal["kOmega"] = "kOmega"
    turbulence: Literal["on", "off"] = "on"
    printCoeffs: Literal["on", "off"] = "on"
    beta: float = 0.075
    gamma: float = 5.0 / 9.0
    sigma_k: float = 2.0
    sigma_omega: float = 2.0


class SpalartAllmaras(IOModelBase):
    RASModel: Literal["SpalartAllmaras"] = "SpalartAllmaras"
    turbulence: Literal["on", "off"] = "on"
    printCoeffs: Literal["on", "off"] = "on"
    Cb1: float = 0.1355
    Cb2: float = 0.622
    sigma: float = 2 / 3
    kappa: float = 0.41

class RASConfig(
    RootModel[
        Annotated[
            Union[kEpsilon, kOmega, SpalartAllmaras],
            Field(discriminator="RASModel"),
        ]
    ],
    IOModelMixin
):
    @classmethod
    def from_ofdict(cls, d):
        """Parse RAS config from OpenFOAM dictionary."""
        ras_model_name = d.get[str]("RASModel")

        # Map model names to classes
        model_classes = {
            "kEpsilon": kEpsilon,
            "kOmega": kOmega, 
            "SpalartAllmaras": SpalartAllmaras
        }
        
        if ras_model_name not in model_classes:
            raise ValueError(f"Unknown RAS model: {ras_model_name}")
        
        model_class = model_classes[ras_model_name]
        model_instance = model_class.from_ofdict(d)
        
        return cls(model_instance)