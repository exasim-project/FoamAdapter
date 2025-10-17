from foamadapter.modules.fields import Fields, Field
import pydantic
from pybFoam import scalarField, vectorField
from foamadapter.modules.setup import visualize_dag

# Factory classes with dependency annotation using Fields.deps
@Fields.deps()
class CreateTemperatureField:
    def __call__(self) -> Field:
        return Field(
            type="scalarField",
            value=scalarField([300, 310, 320]),
            dimensions=(0, 0, 0, 1, 0, 0, 0),
            description="Temperature field"
        )

@Fields.deps("pressure", "temperature")
class CreateDensityField:
    def __call__(self) -> Field:
        return Field(
            type="scalarField",
            value=scalarField([1.225, 1.200, 1.180]),
            dimensions=(1, -3, 0, 0, 0, 0, 0),
            description="Air density field"
        )

@Fields.deps("temperature")
def create_viscosity() -> Field:
    return Field(
        type="scalarField",
        value=scalarField([1.8e-5, 1.8e-5, 1.8e-5]),
        dimensions=(1, -3, 0, 0, 0, 0, 0),
        description="Air viscosity field"
    )

def test_fields():
    fields = Fields()
    assert list(fields.names()) == []

    # Test 1: Add fields before initialization
    velocity_field = Field(
        type="vectorField",
        value=vectorField([(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
        dimensions=(0, 1, -1, 0, 0, 0, 0),
        description="Velocity field"
    )
    pressure_field = Field(
        type="scalarField",
        value=scalarField([101325, 101300]),
        dimensions=(1, -1, -2, 0, 0, 0, 0),
        description="Pressure field"
    )

    fields.add_field("velocity", velocity_field)
    fields.add_field("pressure", pressure_field)

    # Test 2: Add factory functions
    fields.add_field("temperature", CreateTemperatureField())

    fields.add_field("density", CreateDensityField())
    fields.add_field("viscosity", create_viscosity)
    
    # Test 3: Check that fields are not accessible before initialization
    assert list(fields.names()) == ["velocity", "pressure", "temperature", "density", "viscosity"]
    
    visualize_dag(fields.dependencies(), title="Field/Model Dependency DAG", filename=None, show=True)

    # Test 4: Initialize all fields
    fields.initialize_all()
    
    # Test 5: Now fields should be accessible
    assert fields["velocity"].value[0] == [1.0, 0.0, 0.0]
    assert fields["pressure"].value[0] == 101325
    assert fields["temperature"].value[0] == 300
    assert fields["density"].value[0] == 1.225
    assert fields["viscosity"].value[0] == 1.8e-5
    
    # Test 6: Check all fields are resolved
    assert all(isinstance(field, Field) for field in fields.entries.values())
    assert len(fields.entries) == 5  # velocity, pressure, temperature, density, viscosity
    # assert False

    visualize_dag(fields.dependencies(), title="Field/Model Dependency DAG", filename=None, show=True)
