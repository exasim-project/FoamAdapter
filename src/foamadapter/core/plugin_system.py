from pydantic import BaseModel, Field, create_model
from typing import Annotated, Type, Any, List, Dict, Union

extensible_model_registry: Dict[str, Dict] = {}

def register_extensible_model(key: str, base: Type[BaseModel], discriminator: str = "plugin_type"):
    registry: List[Type[BaseModel]] = []
    def register(cls: Type[BaseModel]):
        registry.append(cls)
        extensible_model_registry[key]["model"] = build_model()
        return cls
    def build_model():
        if not registry:
            union = object
        else:
            union = Annotated[Union[tuple(registry)], Field(discriminator=discriminator)]
        return create_model(
            f"{key}Model",
            plugin=(union, ...),
            __base__=base
        )
    extensible_model_registry[key] = {
        "registry": registry,
        "register": register,
        "build_model": build_model,
        "model": build_model()
    }
    return register

def get_extensible_model(key: str):
    return extensible_model_registry[key]["model"]

def register_plugin(key: str, plugin_cls: Type[BaseModel]):
    return extensible_model_registry[key]["register"](plugin_cls)

def list_registered_plugins(key: str):
    return [cls.__name__ for cls in extensible_model_registry[key]["registry"]]

def get_plugin_schema(key: str):
    return get_extensible_model(key).model_json_schema()
