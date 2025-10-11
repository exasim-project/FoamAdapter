from foamadapter.modules.fields import Fields, Field
from pybFoam import scalarField, vectorField


def create_temperature_field():
    """Factory function that creates a temperature field."""
    return Field(
        type="scalarField",
        value=scalarField([300, 310, 320]),
        dimensions=(0, 0, 0, 1, 0, 0, 0),  # Temperature dimensions
        description="Temperature field"
    )


def test_fields():
    fields = Fields({})
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
    fields.add_field("temperature", create_temperature_field)
    
    fields["density"] = lambda: Field(
        type="scalarField",
        value=scalarField([1.225, 1.200, 1.180]),
        dimensions=(1, -3, 0, 0, 0, 0, 0),
        description="Air density field"
    )

    fields["viscosity"] = lambda: Field(
        type="scalarField",
        value=scalarField([1.8e-5, 1.9e-5]),
        dimensions=(1, -1, -1, 0, 0, 0, 0),
        description="Dynamic viscosity"
    )
    
    # Test 3: Check that fields are not accessible before initialization
    assert list(fields.names()) == ["velocity", "pressure", "temperature", "density", "viscosity"]
    
    try:
        _ = fields["velocity"]
        assert False, "Should have raised RuntimeError"
    except RuntimeError as e:
        assert "Fields must be initialized before access" in str(e)
    
    # Test 4: Initialize all fields
    fields.initialize_all()
    
    # Test 5: Now fields should be accessible
    assert fields["velocity"].value[0] == [1.0, 0.0, 0.0]
    assert fields["pressure"].value[0] == 101325
    assert fields["temperature"].value[0] == 300
    assert fields["density"].value[0] == 1.225
    assert fields["viscosity"].value[0] == 1.8e-5
    
    # Test 6: Check all fields are resolved
    assert all(isinstance(field, Field) for field in fields.root.values())
    assert len(fields.root) == 5  # velocity, pressure, temperature, density, viscosity