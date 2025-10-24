
from pydantic import BaseModel, Field, create_model
from typing import Annotated, Type, Union
from dataclasses import dataclass, field


@dataclass
class PluginRegistry:
    base_cls: Type[BaseModel]
    discriminator_variable: str
    discriminator: str
    plugin_registry: list = field(default_factory=list)
    plugin_model: Type[BaseModel] = None

class PluginSystem:
    """
    PluginSystem provides a runtime-extensible plugin/config system using Pydantic discriminated unions.

    Features:
    - Central registry (_registry) for all plugin base types and their plugin classes.
    - Decorator API for explicit registration of plugin families and plugin config classes.
    - Supports multiple independent plugin families, each with its own registry and extensible model.
    - Uses a PluginRegistry dataclass to store metadata for each plugin base type:
        - base_cls: The plugin base class (usually a Pydantic model).
        - plugin_registry: List of registered plugin config classes for this type.
        - discriminator_variable: Name of the field holding the union (e.g., 'plugin').
        - discriminator: Name of the discriminator field in plugin configs (e.g., 'plugin_type').
        - extensible_model: The dynamically generated Pydantic model for this plugin type.

    Usage:
    1. Decorate your base class:
        @PluginSystem.register(discriminator_variable="plugin", discriminator="plugin_type")
        class PluginBase(BaseModel):
            name: str

    2. Register plugin config classes:
        @PluginBase.register
        class AddOneConfig(BaseModel):
            plugin_type: Literal["add_one"]
            amount: int

    3. Access the extensible model:
        PluginModel = PluginBase._extensible_model

    4. The registry can be queried for all plugin families and their plugins:
    PluginSystem._registry["PluginBase"].plugin_registry

    This design allows runtime extensibility, developer-friendly registration, and schema validation for plugin/config systems.
    """
    _registry = {}

    @staticmethod
    def register(discriminator_variable, discriminator):
        def base_decorator(base_cls):
            # Store metadata for this base class in the registry as a dataclass
            PluginSystem._registry[base_cls.__name__] = PluginRegistry(
                base_cls=base_cls,
                plugin_registry=[],
                discriminator_variable=discriminator_variable,
                discriminator=discriminator,
                plugin_model=None,
            )
            # define plugin decorator that is add as register method to base class
            def _register(plugin_cls):
                registry_obj = PluginSystem._registry[base_cls.__name__]
                registry_obj.plugin_registry.append(plugin_cls)
                registry = registry_obj.plugin_registry
                union = Annotated[Union[tuple(registry)], Field(discriminator=registry_obj.discriminator)]
                model = create_model(
                    f"{base_cls.__name__}ExtensibleModel",
                    **{registry_obj.discriminator_variable: (union, ...)},
                    __base__=base_cls
                )
                registry_obj.plugin_model = model
                base_cls.plugin_model = model
                return plugin_cls
            base_cls.register = _register
            # Initial model with no plugins
            registry_obj = PluginSystem._registry[base_cls.__name__]
            union = object
            model = create_model(
                f"{base_cls.__name__}ExtensibleModel",
                **{registry_obj.discriminator_variable: (union, ...)},
                __base__=base_cls
            )
            registry_obj.plugin_model = model
            base_cls.plugin_model = model
            # Add a classmethod 'create' for user-friendly instantiation
            def create(cls, **kwargs):
                return cls.plugin_model(**kwargs)
            base_cls.create = classmethod(create)
            return base_cls
        return base_decorator
