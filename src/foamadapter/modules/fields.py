from pydantic import BaseModel, RootModel, ConfigDict
from typing import Any, Callable, Optional, Union

class Field(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    
    type: str
    value: Any
    dimensions: tuple[int, int, int, int, int, int, int]
    description: str

class Fields(RootModel[dict[str, Union[Field, Callable[[], Field]]]]):
    """Container for OpenFOAM fields that supports both Field objects and factory functions."""
    
    def __init__(self, root=None):
        super().__init__(root or {})
        self._initialized = False
    
    def __getitem__(self, name: str) -> Field:
        if not self._initialized:
            raise RuntimeError("Fields must be initialized before access. Call initialize_all() first.")
        
        item = self.root[name]
        return item
    
    def __setitem__(self, name: str, value: Union[Field, Callable[[], Field]]) -> None:
        """Set a field directly (Field object or factory function)."""
        self.root[name] = value
        
    def add_field(self, name: str, field: Union[Field, Callable[[], Field]]) -> None:
        """Add a field or a function that creates a field."""
        self.root[name] = field
        
    def names(self):
        return self.root.keys()
        
    def keys(self):
        return self.root.keys()
        
    def values(self):
        return self.root.values()
        
    def items(self):
        return self.root.items()
    
    def initialize_all(self) -> None:
        """Resolve all factory functions to Field objects."""
        for name in list(self.root.keys()):  # Use list() to avoid modification during iteration
            item = self.root[name]
            
            # If it's a callable, resolve it
            if callable(item):
                field = item()
                if not isinstance(field, Field):
                    raise ValueError(f"Function must return a Field object, got {type(field)}")
                self.root[name] = field
        
        self._initialized = True
    
