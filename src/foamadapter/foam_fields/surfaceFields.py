import pybFoam


class surfaceScalarField(pybFoam.surfaceScalarField):
    """Wrapper for pybFoam.surfaceScalarField with additional metadata."""

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
        return "surfaceScalarField"


class surfaceVectorField(pybFoam.surfaceVectorField):
    """Wrapper for pybFoam.surfaceVectorField with additional metadata."""

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
        return "surfaceVectorField"
