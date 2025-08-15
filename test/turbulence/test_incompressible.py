from foamadapter.turbulence.incompressible import TurbulenceModel



def test_laminar():
    config = TurbulenceModel.from_file("test/turbulence/incompressible/turbulenceProperties_laminar")
    assert config.simulationType == "laminar"


def test_RAS_komega():
    config = TurbulenceModel.from_file("test/turbulence/incompressible/turbulenceProperties_kOmega")
    assert config.simulationType == "RAS"
    assert config.RAS.root.RASModel == "kOmega"


def test_LES_Smagorinsky():
    config = TurbulenceModel.from_file("test/turbulence/incompressible/turbulenceProperties_Smagorinsky")
    assert config.simulationType == "LES"
    assert config.LES.root.LESModel == "Smagorinsky"
    assert config.LES.root.Cs == 0.17
