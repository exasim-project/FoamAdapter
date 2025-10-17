from pydantic import BaseModel
import pydantic
from typing import Union
from typing import Protocol, runtime_checkable


class Model(BaseModel):
    model_config = {"arbitrary_types_allowed": True}
    type: str
    parameters: dict[str, Union[float, int, str, bool]]
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
    entries: dict[str, Union[Model, ModelsFactory]] = pydantic.Field(default_factory=dict)

    def __init__(self, models: dict):
        self.models = models

    def add_model(self, name: str, model: any) -> None:
        self.models[name] = model

    def __getitem__(self, name: str) -> any:
        return self.models[name]
    

    def deps()