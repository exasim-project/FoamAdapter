from typing import Protocol
from pybFoam import Time, computeCFLNumber, Info



class StabilityCriteria:
    """
    This class is used to check the stability criteria for the pressure-velocity coupling in incompressible flows.
    It ensures that the numerical scheme remains stable during the simulation.
    """

    def __init__(self):
        self.criteria = []
        self.GREAT = 1e30
        self.SMALL = 1e-15
        self.maxDeltaT = 1e-3  # Default maximum deltaT

    def add_criteria(self, criterion):
        self.criteria.append(criterion)

    def setDelta(self, runTime):
        deltaT = runTime.deltaTValue()
        
        if not self.criteria:
            return

        ratios = [criterion.get_stability_ratio() for criterion in self.criteria]

        # limit to 1.2 to avoid too large time steps
        ratios = [min(ratio, 1.2) for ratio in ratios if ratio > self.SMALL and ratio < self.GREAT]

        # Set most restrictive time step
        finalDeltaT = min(min(deltaT * ratio for ratio in ratios), self.maxDeltaT)
        runTime.setDeltaT(finalDeltaT)
        Info(f"deltaT = {runTime.deltaTValue()}")
        runTime.increment()


class StabilityCriterion(Protocol):
    """
    Protocol for stability criteria.
    Each criterion must implement the get_stability_ratio method.
    """

    def get_stability_ratio(self) -> float:
        """
        Get the stability ratio for the criterion.
        :return: Stability ratio as a float.
        """
        pass


class CFLNumber:
    """
    Class to represent the CFL number stability criterion.
    """

    def __init__(self, phi, max_cfl_number: float):
        self.phi = phi
        self.max_cfl_number = max_cfl_number

    def get_stability_ratio(self) -> float:
        """
        Get the stability ratio for the CFL number criterion.
        :return: Stability ratio as a float.
        """
        max_cfl_number, mean_cfl_number = computeCFLNumber(self.phi)
        Info(f"Courant Number mean: {mean_cfl_number}, max: {max_cfl_number}")
        return self.max_cfl_number / max_cfl_number