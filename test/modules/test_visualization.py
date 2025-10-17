"""
Tests for the DAG visualization functions.
"""
import pytest
from pathlib import Path
from foamadapter.modules.fields import Fields, Field, get_field
from foamadapter.modules.models import Models, Model, get_model
from foamadapter.modules.setup import visualize_dag, visualize_containers, initialize_containers
from pybFoam import scalarField, vectorField


# Test factory classes and functions
@Fields.deps()
class CreateTemperatureField:
    def __call__(self, deps: dict = None) -> Field:
        return Field(
            type="scalarField",
            value=scalarField([300, 310, 320]),
            dimensions=(0, 0, 0, 1, 0, 0, 0),
            description="Temperature field"
        )


@Fields.deps("pressure", "temperature")
class CreateDensityField:
    def __call__(self, deps: dict = None) -> Field:
        temperature = get_field(deps, "temperature")
        pressure = get_field(deps, "pressure")
        return Field(
            type="scalarField",
            value=scalarField([1.225, 1.200, 1.180]),
            dimensions=(1, -3, 0, 0, 0, 0, 0),
            description=f"Air density field (T={temperature.description}, P={pressure.description})"
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


def create_test_containers():
    """Create test containers with dependencies for visualization tests."""
    fields = Fields()
    models = Models()
    
    # Add direct field entries
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
    
    # Add field factories
    fields.add_field("temperature", CreateTemperatureField())
    fields.add_field("density", CreateDensityField())
    
    # Add model factories
    models.add_model("viscosity", viscosity_model)
    models.add_model("turbulence", turbulence_model)
    
    return fields, models


def test_visualize_dag_single_container():
    """Test visualizing a single container dependency graph."""
    fields, _ = create_test_containers()
    
    # Test that visualization works without showing (for CI)
    graph = visualize_dag(
        fields.dependencies(), 
        title="Fields Only Dependency Graph",
        show=False,
        filename=None
    )
    
    # Verify graph structure
    assert len(graph.nodes) == 4  # velocity, pressure, temperature, density
    assert len(graph.edges) == 2  # density depends on pressure + temperature
    
    # Verify specific dependencies
    assert "temperature" in graph.nodes
    assert "density" in graph.nodes
    assert ("pressure", "density") in graph.edges
    assert ("temperature", "density") in graph.edges


def test_visualize_dag_basic_styling():
    """Test basic DAG visualization styling options."""
    fields, _ = create_test_containers()
    
    graph = visualize_dag(
        fields.dependencies(),
        title="Fields Basic Styling",
        show=False
    )
    
    assert "density" in graph.nodes
    assert len(graph.nodes) == 4


def test_visualize_containers_multi_container():
    """Test visualizing multiple containers with different colors."""
    fields, models = create_test_containers()
    
    # Test multi-container visualization BEFORE initialization to see dependencies
    graph = visualize_containers(
        fields, models,
        title="Multi-Container Dependency Graph", 
        show=False,
        top_to_bottom=True,
        level_lines=True
    )
    
    # Verify that both container types are represented
    node_names = list(graph.nodes)
    
    # Should have namespaced nodes
    field_nodes = [n for n in node_names if n.startswith("fields:")]
    model_nodes = [n for n in node_names if n.startswith("models:")]
    
    assert len(field_nodes) == 4  # velocity, pressure, temperature, density
    assert len(model_nodes) == 2  # viscosity, turbulence
    
    # Verify cross-container dependencies exist
    edges = list(graph.edges)
    cross_container_edges = [(u, v) for u, v in edges 
                            if u.split(":")[0] != v.split(":")[0]]
    assert len(cross_container_edges) > 0  # Should have field->model dependencies


def test_visualize_dag_save_to_file(tmp_path):
    """Test saving visualization to file."""
    fields, _ = create_test_containers()
    
    output_file = tmp_path / "test_dag.png"
    
    graph = visualize_dag(
        fields.dependencies(),
        title="Test DAG Save",
        filename=output_file,
        show=False
    )
    
    # Verify file was created
    assert output_file.exists()
    assert output_file.stat().st_size > 0  # File has content
    
    # Verify graph is returned
    assert len(graph.nodes) == 4


def test_visualize_containers_save_to_file(tmp_path):
    """Test saving multi-container visualization to file."""
    fields, models = create_test_containers()
    
    output_file = tmp_path / "test_multi_container.png"
    
    graph = visualize_containers(
        fields, models,
        title="Multi-Container Test",
        filename=output_file,
        show=False
    )
    
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    assert len(graph.nodes) == 6  # 4 fields + 2 models


def test_visualize_dag_layout_options():
    """Test different layout and styling options."""
    fields, _ = create_test_containers()
    
    # Test bottom-to-top layout
    graph1 = visualize_dag(
        fields.dependencies(),
        title="Bottom to Top",
        show=False,
        top_to_bottom=False
    )
    
    # Test without level lines
    graph2 = visualize_dag(
        fields.dependencies(),
        title="No Level Lines",
        show=False,
        level_lines=False
    )
    
    # Both should have same structure, just different visual layout
    assert len(graph1.nodes) == len(graph2.nodes) == 4
    assert len(graph1.edges) == len(graph2.edges) == 2


def test_empty_dependencies():
    """Test visualization with empty or minimal dependencies."""
    # Test empty dependencies
    empty_graph = visualize_dag(
        {},
        title="Empty Graph",
        show=False
    )
    assert len(empty_graph.nodes) == 0
    assert len(empty_graph.edges) == 0
    
    # Test single node with no dependencies
    single_node = visualize_dag(
        {"single": set()},
        title="Single Node",
        show=False
    )
    assert len(single_node.nodes) == 1
    assert len(single_node.edges) == 0


if __name__ == "__main__":
    # Manual test with visualization display (for development)
    print("Creating test containers...")
    fields, models = create_test_containers()
    
    print("\n1. Single container (Fields only):")
    visualize_dag(
        fields.dependencies(),
        title="Fields Dependency Graph",
        show=True
    )
    
    print("\n2. Multi-container visualization:")
    initialize_containers(fields, models)
    
    visualize_containers(
        fields, models,
        title="Multi-Container Dependencies",
        show=True,
        top_to_bottom=True,
        level_lines=True
    )
    
    print("Manual tests completed!")