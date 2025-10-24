from typing import Literal, Union, Annotated, Optional
from pydantic import Field, RootModel, create_model
from pybFoam.io.model_base import IOModelBase, IOModelMixin
from foamadapter.inputs_files.case_inputs import Registry, FileSpec


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

    @classmethod
    def additional_inputs(cls, reg: Registry) -> Registry:
        """Example method for kEpsilon properties."""
        FvSchemes, fileSpec = reg["fvSchemes"]
        divSchemes = FvSchemes.model_fields["divSchemes"].annotation
        updated_divSchemes = create_model(
            "divSchemes",
            __base__=divSchemes,
            div_phi_k=(Optional[str], Field(alias="div(phi,k)", default=None)),
            div_phi_epsilon=(str, Field(alias="div(phi,epsilon)")),
        )
        
        update_fvSchemes = create_model(
            "FvSchemes",
            __base__=FvSchemes,
            divSchemes=(updated_divSchemes, ...) 
        )

        reg["fvSchemes"] = (update_fvSchemes, fileSpec)
        return reg


class kOmega(IOModelBase):
    RASModel: Literal["kOmega"] = "kOmega"
    turbulence: Literal["on", "off"] = "on"
    printCoeffs: Literal["on", "off"] = "on"
    beta: float = 0.075
    gamma: float = 5.0 / 9.0
    sigma_k: float = 2.0
    sigma_omega: float = 2.0

    @classmethod
    def additional_inputs(cls, reg: Registry) -> Registry:
        """Example method for kOmega properties."""

        FvSchemes, fileSpec = reg["fvSchemes"]
        divSchemes = FvSchemes.model_fields["divSchemes"].annotation
        updated_divSchemes = create_model(
            "divSchemes",
            __base__=divSchemes,
            div_phi_k=(Optional[str], Field(alias="div(phi,k)", default=None)),
            div_phi_omega=(str, Field(alias="div(phi,omega)")),
        )
        
        update_fvSchemes = create_model(
            "FvSchemes",
            __base__=FvSchemes,
            divSchemes=(updated_divSchemes, ...) 
        )

        reg["fvSchemes"] = (update_fvSchemes, fileSpec)

        return reg


class SpalartAllmaras(IOModelBase):
    RASModel: Literal["SpalartAllmaras"] = "SpalartAllmaras"
    turbulence: Literal["on", "off"] = "on"
    printCoeffs: Literal["on", "off"] = "on"
    Cb1: float = 0.1355
    Cb2: float = 0.622
    sigma: float = 2 / 3
    kappa: float = 0.41

    @classmethod
    def additional_inputs(cls, reg: Registry) -> Registry:
        """Example method for SpalartAllmaras properties."""
        # reg["SpalartAllmaras"] = 2
        return reg

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
    
    @classmethod
    def additional_inputs(cls, d, reg: Registry) -> Registry:
        """Example method for RAS properties."""

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

        return model_class.additional_inputs(reg)