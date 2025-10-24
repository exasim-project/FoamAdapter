from foamadapter.inputs_files.case_inputs import Registry, FileSpec
from foamadapter.turbulence.incompressible import TurbulenceModel
from foamadapter.inputs_files.system import FvSchemesBase


def test_laminar():
    test_file = "test/turbulence/incompressible/turbulenceProperties_laminar"
    config = TurbulenceModel.from_file(test_file)
    assert config.simulationType == "laminar"

    reg = {

    }
    reg = TurbulenceModel.inputs(reg, turb_file=test_file)
    assert reg["turbulenceProperties"] == (
        TurbulenceModel,
        FileSpec(test_file, encoding="utf-8", required=True)
    )


def test_RAS_komega():
    test_file = "test/turbulence/incompressible/turbulenceProperties_kOmega"
    config = TurbulenceModel.from_file(test_file)
    assert config.simulationType == "RAS"
    assert config.RAS.root.RASModel == "kOmega"

    reg = {
        "fvSchemes": (FvSchemesBase, FileSpec("system/fvSchemes", encoding="utf-8", required=True))
    }
    reg = TurbulenceModel.inputs(reg, turb_file=test_file)
    assert reg["turbulenceProperties"] == (
        TurbulenceModel,
        FileSpec(test_file, encoding="utf-8", required=True)
    )
    # Direct access to divSchemes field
    FvSchemes = reg["fvSchemes"][0]
    divSchemes_field = FvSchemes.model_fields["divSchemes"]
    divSchemes_type = divSchemes_field.annotation    
   
    # Check if specific fields exist and get their aliases
    assert "div_phi_omega" in divSchemes_type.model_fields
    assert "div_phi_k" in divSchemes_type.model_fields
    
    assert divSchemes_type.model_fields["div_phi_omega"].alias == "div(phi,omega)"
    assert divSchemes_type.model_fields["div_phi_k"].alias == "div(phi,k)"



def test_LES_Smagorinsky():
    test_file = "test/turbulence/incompressible/turbulenceProperties_Smagorinsky"
    config = TurbulenceModel.from_file(test_file)
    assert config.simulationType == "LES"
    assert config.LES.root.LESModel == "Smagorinsky"
    assert config.LES.root.Cs == 0.17

    reg = {
    }
    reg = TurbulenceModel.inputs(reg, turb_file=test_file)
    assert reg["turbulenceProperties"] == (
        TurbulenceModel,
        FileSpec(test_file, encoding="utf-8", required=True)
    )