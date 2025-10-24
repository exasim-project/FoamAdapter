from pydantic import BaseModel, RootModel, ConfigDict
import pydantic
from typing import Any, Callable, Optional, Union, Type, Generic, TypeVar
from typing import Protocol, runtime_checkable


@runtime_checkable
class Field(Protocol):

    @property
    def type(self) -> str: ...

    @property
    def dimensions(self) -> tuple[int, int, int, int, int, int, int]: ...

    @property
    def description(self) -> str: ...


@runtime_checkable
class FieldFactory(Protocol):
    """Callable that creates a Field object."""

    @property
    def dependencies(self) -> list[str]: ...

    def __call__(self, deps: dict[str, Any]) -> Field: ...


# class Fields(RootModel[dict[str, Union[Field, FieldFactory]]]):
class Fields(BaseModel):
    @staticmethod
    def deps(*dependencies: str):
        """Decorator to attach dependencies to a FieldFactory, setting .dependencies for compatibility."""

        def decorator(factory):
            deps_list = list(dependencies)
            setattr(factory, "dependencies", deps_list)
            return factory

        return decorator

    """Container for OpenFOAM fields that supports both Field objects and factory functions."""
    model_config = {"arbitrary_types_allowed": True}
    entries: dict[str, Union[Field, FieldFactory]] = pydantic.Field(
        default_factory=dict
    )

    def __getitem__(self, name: str) -> Field:
        item = self.entries[name]
        return item

    def __setitem__(self, name: str, value: Union[Field, FieldFactory]) -> None:
        """Set a field directly (Field object or factory function)."""
        # Validate that the field implements the correct protocol
        if not isinstance(value, Field) and not isinstance(value, FieldFactory):
            raise TypeError(
                f"Expected Field or FieldFactory for '{name}', got {type(value).__name__}"
            )
        self.entries[name] = value

    def add_field(self, name: str, field: Union[Field, FieldFactory]) -> None:
        """Add a field or a function that creates a field."""
        # Validate that the field implements the correct protocol
        if not isinstance(field, Field) and not isinstance(field, FieldFactory):
            raise TypeError(
                f"Expected Field or FieldFactory for '{name}', got {type(field).__name__}"
            )
        self.entries[name] = field

    def add_fields(self, fields: dict[str, Union[Field, FieldFactory]]) -> None:
        """Add multiple fields from a dictionary."""
        for name, field in fields.items():
            self.add_field(name, field)  # Use add_field to get validation

    def dependencies(self) -> dict[str, set[str]]:
        deps = {}
        for name, item in self.entries.items():
            if callable(item) and hasattr(item, "dependencies"):
                deps[name] = set(item.dependencies)
            else:
                deps[name] = set()
        return deps

    def names(self):
        return self.entries.keys()

    def keys(self):
        return self.entries.keys()

    def values(self):
        return self.entries.values()

    def items(self):
        return self.entries.items()

    def is_initialized(self) -> bool:
        return all(not callable(item) for item in self.entries.values())


from ..foam_fields.volFields import volScalarField, volVectorField

volFields = TypeVar("T", volScalarField, volVectorField)


class FieldReader(Generic[volFields]):
    """Generic field reader that works as a callable with Fields."""

    model_config = {"arbitrary_types_allowed": True}

    def __init__(
        self,
        field_class: Type[volFields],
        mesh,
        field_name: str,
        description: str,
        dimensions: tuple[int, int, int, int, int, int, int],
    ):
        self.field_class = field_class
        self.mesh = mesh
        self.field_name = field_name
        self._description = description
        self._dimensions = dimensions

    @property
    def dependencies(self) -> list[str]:
        return []

    def __call__(self, deps: dict) -> Field:
        """Read field from disk and convert to Field object."""
        try:
            foam_field = self.field_class.read_field(self.mesh, self.field_name)
            return self.field_class(
                value=foam_field,
                dimensions=self._dimensions,  # Placeholder; implement dimension extraction if needed
                description=self._description,
            )
        except Exception as e:
            raise RuntimeError(
                f"Failed to read {self.field_class.__name__} '{self.field_name}': {e}"
            )


def get_field(deps: dict, key: str) -> Field:
    """Retrieve a Field dependency with type validation."""
    obj = deps[key]
    if not isinstance(obj, Field):
        raise TypeError(
            f"Expected Field for dependency '{key}', got {type(obj).__name__}"
        )
    return obj
