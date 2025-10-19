"""
Test for an extensible plugin config system using pydantic discriminated unions and a registry pattern.
Refactored to use a generic registry and factory for multiple extensible models.
"""

from pydantic import BaseModel, Field, ValidationError
from typing import Literal
import pytest
from foamadapter.core.plugin_system import register_extensible_model, get_extensible_model, register_plugin

# --- Example usage ---
class PluginBase(BaseModel):
    name: str

# Get the decorator for plugin registration
plugin_decorator = register_extensible_model("Plugin", PluginBase, discriminator="plugin_type")

@plugin_decorator
class AddOneConfig(BaseModel):
    plugin_type: Literal["add_one"]
    amount: int

@plugin_decorator
class SquareConfig(BaseModel):
    plugin_type: Literal["square"]
    factor: int

@plugin_decorator
class DoubleConfig(BaseModel):
    plugin_type: Literal["double"]
    multiplier: int

# Extensible types (added in tests)
class TripleConfig(BaseModel):
    plugin_type: Literal["triple"]
    factor: int

class HamsterConfig(BaseModel):
    plugin_type: Literal["hamster"]
    runs_in_wheel: bool

PluginModel = get_extensible_model("Plugin")

def test_valid_configs_n_models():
    m1 = PluginModel(plugin={"plugin_type": "add_one", "amount": 2}, name="adder")
    assert m1.plugin.amount == 2
    m2 = PluginModel(plugin={"plugin_type": "square", "factor": 3}, name="squarer")
    assert m2.plugin.factor == 3
    m3 = PluginModel(plugin={"plugin_type": "double", "multiplier": 4}, name="doubler")
    assert m3.plugin.multiplier == 4

def test_invalid_config_n_models():
    with pytest.raises(ValidationError):
        PluginModel(plugin={"plugin_type": "add_one"}, name="bad")  # missing 'amount'
    with pytest.raises(ValidationError):
        PluginModel(plugin={"plugin_type": "square", "amount": 2}, name="bad")  # wrong field

def test_extensibility_n_models():
    register_plugin("Plugin", TripleConfig)
    PluginModel = get_extensible_model("Plugin")
    m4 = PluginModel(plugin={"plugin_type": "triple", "factor": 5}, name="tripler")
    assert m4.plugin.factor == 5
    assert m4.plugin.plugin_type == "triple"

    register_plugin("Plugin", HamsterConfig)
    PluginModel = get_extensible_model("Plugin")
    m5 = PluginModel(plugin={"plugin_type": "hamster", "runs_in_wheel": True}, name="hammy")
    assert m5.plugin.runs_in_wheel is True
    assert m5.plugin.plugin_type == "hamster"

def test_json_schema():
    PluginModel = get_extensible_model("Plugin")
    schema = PluginModel.model_json_schema()
    discriminator = schema["properties"]["plugin"]["discriminator"]
    assert discriminator["propertyName"] == "plugin_type"
    mapping = discriminator["mapping"]
    assert "add_one" in mapping
    assert "square" in mapping
    assert "double" in mapping
    register_plugin("Plugin", TripleConfig)
    PluginModel = get_extensible_model("Plugin")
    schema = PluginModel.model_json_schema()
    mapping = schema["properties"]["plugin"]["discriminator"]["mapping"]
    assert "triple" in mapping
    register_plugin("Plugin", HamsterConfig)
    PluginModel = get_extensible_model("Plugin")
    schema = PluginModel.model_json_schema()
    mapping = schema["properties"]["plugin"]["discriminator"]["mapping"]
    assert "hamster" in mapping
