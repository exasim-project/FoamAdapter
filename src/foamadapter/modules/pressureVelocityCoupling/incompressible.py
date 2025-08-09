from pybFoam import (
    volScalarField, volVectorField, surfaceScalarField, fvScalarMatrix, fvVectorMatrix,
    fvMesh, Time, fvc, fvm, Word, dictionary, Info,
    solve, adjustPhi, constrainPressure, createPhi, setRefCell,
    constrainHbyA, pimpleControl, computeCFLNumber
)
from ..stability_criteria import CFLNumber


class PimpleAlgorithm:
    """
    Protocol for PIMPLE algorithm implementations.
    This can be used to define the interface for different PIMPLE algorithms.
    """

    def __init__(self, mesh, turbulence=None, laminarTransport=None):
        """
        Initialize the PIMPLE algorithm with the mesh and optional turbulence model.
        :param mesh: The computational mesh.
        :param turbulence: Optional turbulence model.
        :param laminarTransport: Optional laminar transport model.
        """
        self.mesh = mesh
        self.turbulence = turbulence
        self.laminarTransport = laminarTransport
        self.p = volScalarField.read_field(mesh, "p")
        self.U = volVectorField.read_field(mesh, "U")
        self.phi = createPhi(self.U)
        self.UEqn = None  # Placeholder for momentum equation

        fvSolution = dictionary.read("system/fvSolution")
        self.pRefCell, self.pRefValue = setRefCell(self.p, fvSolution.subDict("PIMPLE"))
        self.mesh.setFluxRequired(Word("p"))
        self.cfl_number = CFLNumber(self.phi, 5.0)

    def stability_criteria(self):
        return self.cfl_number

    def momentum_equation(self,pimple) -> None:
        """
        Solve the momentum equations using the PIMPLE algorithm.
        """
        self.UEqn = fvVectorMatrix(
            fvm.ddt(self.U) + fvm.div(self.phi, self.U) + self.turbulence.divDevReff(self.U)
        )

        self.UEqn.relax()

        if pimple.momentumPredictor():
            solve(self.UEqn + fvc.grad(self.p))

    def pressure_correction(self,pimple) -> None:
        """
        Correct the solution based on the PIMPLE algorithm.
        """
        rAU = volScalarField(Word("rAU"), 1.0 / self.UEqn.A())
        HbyA = volVectorField(constrainHbyA(rAU * self.UEqn.H(), self.U, self.p))

        phiHbyA = surfaceScalarField(
            Word("phiHbyA"),
            fvc.flux(HbyA) + fvc.interpolate(rAU) * fvc.ddtCorr(self.U, self.phi)
        )

        adjustPhi(phiHbyA, self.U, self.p)
        constrainPressure(self.p, self.U, phiHbyA, rAU)

        while pimple.correctNonOrthogonal():
            pEqn = fvScalarMatrix(fvm.laplacian(rAU, self.p) - fvc.div(phiHbyA))
            pEqn.setReference(self.pRefCell, self.pRefValue, False)
            pEqn.solve(self.p.select(pimple.finalInnerIter()))

            if pimple.finalNonOrthogonalIter():
                self.phi.assign(phiHbyA - pEqn.flux())

        # Optionally include continuityErrs()
        self.U.assign(HbyA - rAU * fvc.grad(self.p))
        self.U.correctBoundaryConditions()

