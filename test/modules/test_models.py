from foamadapter.modules.models import Models, Model


class DummyEnv:
    def __init__(self, values):
        self._values = values

    def field(self, name):
        return self._values[name]

# Class-based factory
@Models.deps("p0", "T0")
class ViscosityModel:
    def __call__(self, deps):
        p0 = deps["p0"]
        T0 = deps["T0"]
        return f"viscosity({p0}, {T0})"

@Models.deps("U0", "nu")
def turbulenceModel(deps):
    U0 = deps["U0"]
    nu = deps["nu"]
    return f"TurbulenceModel(U0={U0}, nu={nu})"


def test_models_deps_and_registration():
    models = Models()


    models.add_model("turbulenceModel", turbulenceModel)

    # Check dependencies
    assert models["turbulenceModel"].dependencies == ["U0", "nu"]

    # Simulate dependency injection
    env = DummyEnv({"U0": 1.0, "nu": 0.01})
    result = models["turbulenceModel"]({"U0": env.field("U0"), "nu": env.field("nu")})
    assert result == "TurbulenceModel(U0=1.0, nu=0.01)"



    models.add_model("viscosity", ViscosityModel())
    assert models["viscosity"].dependencies == ["p0", "T0"]
    result = models["viscosity"]({"p0": 101325, "T0": 300})
    assert result == "viscosity(101325, 300)"
