"""
Utility field classes for testing that implement the Field protocol.
These classes extend pybFoam field types with the required properties.
"""

from pybFoam import scalarField, vectorField


class RegisteredScalarField(scalarField):
    """
    ScalarField that implements the Field protocol by adding metadata properties.
    Extends pybFoam scalarField with description, dimensions, and type properties.
    """

    def __init__(self, values, description: str, dimensions: tuple[int, int, int, int, int, int, int]):
        super().__init__(values)
        self._description = description
        self._dimensions = dimensions

    @property
    def description(self) -> str:
        return self._description

    @property
    def dimensions(self) -> tuple[int, int, int, int, int, int, int]:
        return self._dimensions

    @property
    def type(self) -> str:
        return "scalarField"


class RegisteredVectorField(vectorField):
    """
    VectorField that implements the Field protocol by adding metadata properties.
    Extends pybFoam vectorField with description, dimensions, and type properties.
    """

    def __init__(self, values, description: str, dimensions: tuple[int, int, int, int, int, int, int]):
        super().__init__(values)
        self._description = description
        self._dimensions = dimensions

    @property
    def description(self) -> str:
        return self._description

    @property
    def dimensions(self) -> tuple[int, int, int, int, int, int, int]:
        return self._dimensions

    @property
    def type(self) -> str:
        return "vectorField"