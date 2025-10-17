from pydantic import BaseModel, RootModel, ConfigDict
import pydantic
from typing import Any, Callable, Optional, Union , Type, Generic, TypeVar
from pybFoam import volScalarField, volVectorField
from typing import Protocol, runtime_checkable

class Field(BaseModel):
    model_config = {"arbitrary_types_allowed": True}
    
    type: str
    value: Any
    dimensions: tuple[int, int, int, int, int, int, int]
    description: str

@runtime_checkable
class FieldFactory(Protocol):
    """Callable that creates a Field object."""
    
    @property
    def dependencies(self) -> list[str]:
        ...

    def __call__(self, deps: dict[str, Any]) -> Field:
        ...

# class Fields(RootModel[dict[str, Union[Field, FieldFactory]]]):
class Fields(BaseModel):
    """Container for OpenFOAM fields that supports both Field objects and factory functions."""
    model_config = {"arbitrary_types_allowed": True}
    entries: dict[str, Union[Field, FieldFactory]] = pydantic.Field(default_factory=dict)
    _initialized: bool = False

    def __getitem__(self, name: str) -> Field:
        if not self._initialized:
            raise RuntimeError("Fields must be initialized before access. Call initialize_all() first.")
        
        item = self.entries[name]
        return item
    
    def __setitem__(self, name: str, value: Union[Field, Callable[[], Field]]) -> None:
        """Set a field directly (Field object or factory function)."""
        self.entries[name] = value
        
    def add_field(self, name: str, field: Union[Field, Callable[[], Field]]) -> None:
        """Add a field or a function that creates a field."""
        self.entries[name] = field

    def add_fields(self, fields: dict[str, Union[Field, Callable[[], Field]]]) -> None:
        """Add multiple fields from a dictionary."""
        for name, field in fields.items():
            self.entries[name] = field

    def dependencies(self) -> dict[str, set[str]]:
        deps = {}
        for name, item in self.entries.items():
            if isinstance(item, FieldFactory):
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
    
    def initialize_all(self) -> None:
        """Resolve all factory functions to Field objects."""
        for name in list(self.entries.keys()):  # Use list() to avoid modification during iteration
            item = self.entries[name]

            # If it's a factory, resolve it
            if isinstance(item, FieldFactory):
                print(f"Initializing field '{name}' using factory.")
                print(f"Dependencies: {item.dependencies}")
                field = item()
                if not isinstance(field, Field):
                    raise ValueError(f"Function must return a Field object, got {type(field)}")
                self.entries[name] = field
        
        self._initialized = True
    
volFields = TypeVar('T', volScalarField, volVectorField)

class FieldReader(Generic[volFields]):
    """Generic field reader that works as a callable with Fields."""
    model_config = {"arbitrary_types_allowed": True}
    
    def __init__(self, field_class: Type[volFields], mesh, field_name: str, dependencies: Optional[list[str]] = None):
        self.field_class = field_class
        self.mesh = mesh
        self.field_name = field_name
        self._dependencies = dependencies or []

    @property
    def dependencies(self) -> list[str]:
        return self._dependencies


    def __call__(self) -> Field:
        """Read field from disk and convert to Field object."""
        try:
            foam_field = self.field_class.read_field(self.mesh, self.field_name)
            return Field(
                type=self.field_class.__name__,
                value=foam_field,
                # dimensions=self._get_dimensions(foam_field),
                dimensions=(0, 0, 0, 0, 0, 0, 0),  # Placeholder; implement dimension extraction if needed
                description=f"{self.field_name} ({self.field_class.__name__})"
            )
        except Exception as e:
            raise RuntimeError(f"Failed to read {self.field_class.__name__} '{self.field_name}': {e}")
    