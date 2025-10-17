from foamadapter.modules.models import Models, Model, get_model
from foamadapter.modules.fields import Fields, Field, get_field
from foamadapter.modules.setup import initialize_containers
from pybFoam import scalarField, vectorField

@Fields.deps()
class CreateTemperatureField:
    def __call__(self, deps: dict ) -> Field:
        return Field(
            type="scalarField",
            value=scalarField([300, 310, 320]),
            dimensions=(0, 0, 0, 1, 0, 0, 0),
            description="Temperature field"
        )
    
@Models.deps("temperature")
def viscosity_model(deps: dict) -> Model:
    temperature = get_field(deps, "temperature")
    return Model(
        type="viscosity",
        parameters={"temperature": temperature},
        description="Viscosity model based on temperature"
    )

@Models.deps("viscosity", "velocity")
def turbulence_model(deps: dict) -> Model:
    viscosity = get_model(deps, "viscosity")
    velocity = get_field(deps, "velocity")
    return Model(
        type="turbulence",
        parameters={"viscosity": viscosity, "velocity": velocity},
        description="Turbulence model based on velocity and viscosity"
    )

def test_initialize_containers():
    fields = Fields()
    models = Models()
    # Add direct entries
    fields.add_field("temperature", CreateTemperatureField())

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

def test_type_checking_helpers():
    """Test that helper functions provide proper type checking."""
    
    # Test get_field with correct type
    field = Field(type="test", value=None, dimensions=(0,0,0,0,0,0,0), description="Test field")
    deps_with_field = {"my_field": field}
    retrieved_field = get_field(deps_with_field, "my_field")
    assert retrieved_field is field
    
    # Test get_model with correct type
    model = Model(type="test", parameters={}, description="Test model")
    deps_with_model = {"my_model": model}
    retrieved_model = get_model(deps_with_model, "my_model")
    assert retrieved_model is model
    
    # Test get_field fails with wrong type
    deps_with_wrong_type = {"my_field": model}  # Model instead of Field
    try:
        get_field(deps_with_wrong_type, "my_field")
        assert False, "Should have raised TypeError"
    except TypeError as e:
        assert "Expected Field for dependency 'my_field', got Model" in str(e)
    
    # Test get_model fails with wrong type
    deps_with_wrong_type = {"my_model": field}  # Field instead of Model
    try:
        get_model(deps_with_wrong_type, "my_model")
        assert False, "Should have raised TypeError"
    except TypeError as e:
        assert "Expected Model for dependency 'my_model', got Field" in str(e)

