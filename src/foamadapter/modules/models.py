from pydantic import BaseModel
import pydantic
from typing import Protocol, Union, Callable, runtime_checkable, Any


class Model(BaseModel):
    model_config = {"arbitrary_types_allowed": True}
    type: str
    parameters: dict[str, Any]
    description: str

@runtime_checkable
class ModelsFactory(Protocol):
    """Callable that creates a Model object."""

    @property
    def dependencies(self) -> list[str]:
        ...

    def __call__(self, deps: dict[str, any]) -> any:
        ...


class Models(BaseModel):
    model_config = {"arbitrary_types_allowed": True}
    entries: dict[str, Union[Model, ModelsFactory]] = pydantic.Field(default_factory=dict)

    def add_model(self, name: str, model: Union[Model, ModelsFactory]) -> None:
        self.entries[name] = model

    def __getitem__(self, name: str) -> Union[Model, ModelsFactory]:
        return self.entries[name]

    def dependencies(self) -> dict[str, set[str]]:
        deps = {}
        for name, item in self.entries.items():
            if callable(item) and hasattr(item, "dependencies"):
                deps[name] = set(item.dependencies)
            else:
                deps[name] = set()
        return deps

    @staticmethod
    def deps(*dependencies: str):
        """Decorator to attach dependencies to a ModelsFactory, setting .dependencies for compatibility."""
        def decorator(factory):
            deps_list = list(dependencies)
            setattr(factory, "dependencies", deps_list)
            return factory
        return decorator


def get_model(deps: dict, key: str) -> Model:
    """Retrieve a Model dependency with type validation."""
    obj = deps[key]
    if not isinstance(obj, Model):
        raise TypeError(f"Expected Model for dependency '{key}', got {type(obj).__name__}")
    return obj