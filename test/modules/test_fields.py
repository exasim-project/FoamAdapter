from foamadapter.modules.fields import Fields, Field, get_field
from foamadapter.modules.setup import initialize_containers
import pydantic
from pybFoam import scalarField, vectorField
from test.modules.utils_fields import RegisteredScalarField, RegisteredVectorField
    
@Fields.deps()
class CreateTemperatureField:
    def __call__(self, deps) -> Field:
        return RegisteredScalarField(
            values=scalarField([300, 310, 320]),
            dimensions=(0, 0, 0, 1, 0, 0, 0),
            description="Temperature field"
        )

@Fields.deps("pressure", "temperature")
class CreateDensityField:
    def __call__(self, deps) -> Field:
        pressure = get_field(deps, "pressure")
        temperature = get_field(deps, "temperature")
        return RegisteredScalarField(
            values=scalarField([1.225, 1.200, 1.180]),
            dimensions=(1, -3, 0, 0, 0, 0, 0),
            description="Air density field"
        )

@Fields.deps("temperature")
def create_viscosity(deps: dict) -> Field:
    temperature = get_field(deps, "temperature")
    return RegisteredScalarField(
        values=scalarField([1.8e-5, 1.8e-5, 1.8e-5]),
        dimensions=(1, -3, 0, 0, 0, 0, 0),
        description=f"Air viscosity field (T={temperature.description})"
    )

def test_fields():
    fields = Fields()
    assert list(fields.names()) == []

    # Add fields before initialization
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

    # Add factory functions
    fields.add_field("temperature", CreateTemperatureField())

    fields.add_field("density", CreateDensityField())
    fields.add_field("viscosity", create_viscosity)
    
    # Check that fields are not accessible before initialization
    assert list(fields.names()) == ["velocity", "pressure", "temperature", "density", "viscosity"]
    
    assert fields.is_initialized() == False  # Before initialization, should be False
    # Initialize all fields
    initialize_containers(fields)

    assert fields.is_initialized() == True
    
    # Now fields should be accessible
    assert fields["velocity"][0] == [1.0, 0.0, 0.0]
    assert fields["pressure"][0] == 101325
    assert fields["temperature"][0] == 300
    assert fields["density"][0] == 1.225
    assert fields["viscosity"][0] == 1.8e-5
    
    # Check all fields are resolved
    assert all(isinstance(field, Field) for field in fields.entries.values())
    assert len(fields.entries) == 5  # velocity, pressure, temperature, density, viscosity
