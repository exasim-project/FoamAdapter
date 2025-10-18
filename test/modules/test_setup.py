from foamadapter.modules.models import Models, Model, get_model
from foamadapter.modules.fields import Fields, Field, get_field
from foamadapter.modules.setup import initialize_containers
from pybFoam import scalarField, vectorField
from test.modules.utils_fields import RegisteredScalarField, RegisteredVectorField

@Fields.deps()
class CreateTemperatureField:
    def __call__(self, deps: dict ) -> Field:
        return RegisteredScalarField(
            values=scalarField([300, 310, 320]),
            dimensions=(0, 0, 0, 1, 0, 0, 0),
            description="Temperature field"
        )
    
@Models.deps("temperature")
def viscosity_model(deps: dict) -> Model:
    temperature = get_field(deps, "temperature")
    class ViscosityModel:

        @property
        def description(self) -> str:
            return "Viscosity model based on temperature"
    return ViscosityModel()

@Models.deps("viscosity", "velocity")
def turbulence_model(deps: dict) -> Model:
    viscosity = get_model(deps, "viscosity")
    velocity = get_field(deps, "velocity")
    class TurbulenceModel:

        @property
        def description(self) -> str:
            return "Turbulence model based on velocity and viscosity"

    return TurbulenceModel()

def test_initialize_containers():
    fields = Fields()
    models = Models()
    # Add direct entries
    fields.add_field("temperature", CreateTemperatureField())

    velocity_field = RegisteredVectorField(
        values=vectorField([(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
        dimensions=(0, 1, -1, 0, 0, 0, 0),
        description="Velocity field"
    )
    pressure_field = RegisteredScalarField(
        values=scalarField([101325, 101300]),
        dimensions=(1, -1, -2, 0, 0, 0, 0),
        description="Pressure field"
    )

    fields.add_field("velocity", velocity_field)
    fields.add_field("pressure", pressure_field)

    models.add_model("viscosity", viscosity_model)
    models.add_model("turbulence", turbulence_model)
    # Add factories with dependencies
    initialize_containers(fields, models)
    # Check that factories are resolved
    assert isinstance(fields.entries["temperature"], Field)
    assert isinstance(models.entries["viscosity"], Model)
    assert isinstance(models.entries["turbulence"], Model)
    assert fields.entries["temperature"].description == "Temperature field"
    assert models.entries["viscosity"].description == "Viscosity model based on temperature"


