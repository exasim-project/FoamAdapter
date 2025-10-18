import pybFoam


class volScalarField(pybFoam.volScalarField):
    """Wrapper for pybFoam.volScalarField with additional metadata."""

    def __init__(
        self,
        value,
        description: str = "",
        dimensions: tuple[int, int, int, int, int, int, int] = (0, 0, 0, 0, 0, 0, 0),
    ):
        super().__init__(value)
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
        return "volScalarField"


class volVectorField(pybFoam.volVectorField):
    """Wrapper for pybFoam.volVectorField with additional metadata."""

    def __init__(
        self,
        value,
        description: str = "",
        dimensions: tuple[int, int, int, int, int, int, int] = (0, 0, 0, 0, 0, 0, 0),
    ):
        super().__init__(value)
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
        return "volVectorField"
