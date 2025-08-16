from foamadapter.inputs_files.case_inputs import Registry, FileSpec
from foamadapter.turbulence.incompressible import TurbulenceModel


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

    }
    reg = TurbulenceModel.inputs(reg, turb_file=test_file)
    assert reg["turbulenceProperties"] == (
        TurbulenceModel,
        FileSpec(test_file, encoding="utf-8", required=True)
    )
    assert reg["kOmega"] == 1


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
    assert reg["Smagorinsky"] == 1