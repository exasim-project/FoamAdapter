from typing import Literal, Optional
from pybFoam import (
    volScalarField,
    volVectorField,
    surfaceScalarField,
    fvScalarMatrix,
    fvVectorMatrix,
    fvMesh,
    Time,
    fvc,
    fvm,
    Word,
    dictionary,
    Info,
    solve,
    adjustPhi,
    constrainPressure,
    createPhi,
    setRefCell,
    constrainHbyA,
    pimpleControl,
    computeCFLNumber,
)
from ..stability_criteria import CFLNumber
from ...inputs_files.case_inputs import Registry
from ...inputs_files.system import ControlDictBase, FvSchemesBase, DIVSchemes
from pydantic import create_model, Field
from ...inputs_files.case_inputs import Registry, FileSpec

ControlDict = create_model(
    "controlDict",
    maxCo=(Optional[float], Field(default=None, description="Maximum Courant number")),
    adjustTimeStep=(Literal["yes", "no"], Field(default="no", description="Adjust time step")),
    __base__=ControlDictBase,
)

divSchemes = create_model(
    "divSchemes",
    __base__=DIVSchemes,
)

FvSchemes = create_model("fvSchemes", __base__=FvSchemesBase)


class PimpleAlgorithm:
    """
    Protocol for PIMPLE algorithm implementations.
    This can be used to define the interface for different PIMPLE algorithms.
    """

    @staticmethod
    def inputs(registry: Registry) -> Registry:
        registry["controlDict"] = (
            ControlDict,
            FileSpec("system/controlDict", description="Time control settings"),
        )
        registry["fvSchemes"] = (
            FvSchemes,
            FileSpec("system/fvSchemes", description="Discretization schemes"),
        )

        return registry

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

    def stability_criteria(self) -> CFLNumber:
        return self.cfl_number

    def momentum_equation(self, pimple) -> None:
        """
        Solve the momentum equations using the PIMPLE algorithm.
        """
        self.UEqn = fvVectorMatrix(
            fvm.ddt(self.U)
            + fvm.div(self.phi, self.U)
            + self.turbulence.divDevReff(self.U)
        )

        self.UEqn.relax()

        if pimple.momentumPredictor():
            solve(self.UEqn + fvc.grad(self.p))

    def pressure_correction(self, pimple) -> None:
        """
        Correct the solution based on the PIMPLE algorithm.
        """
        rAU = volScalarField(Word("rAU"), 1.0 / self.UEqn.A())
        HbyA = volVectorField(constrainHbyA(rAU * self.UEqn.H(), self.U, self.p))

        phiHbyA = surfaceScalarField(
            Word("phiHbyA"),
            fvc.flux(HbyA) + fvc.interpolate(rAU) * fvc.ddtCorr(self.U, self.phi),
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


    # def pressure_correction(self, pimple) -> None:
    #     """
    #     Correct the solution based on the PIMPLE algorithm (with buoyancy, MRF, and non-orthogonal corrections).
    #     """
    #     # Reciprocal of UEqn diagonal
    #     rAU = volScalarField(Word("rAU"), 1.0 / self.UEqn.A())
    #     rAUf = surfaceScalarField(Word("rAUf"), fvc.interpolate(rAU))

    #     # Momentum predictor HbyA
    #     HbyA = volVectorField(constrainHbyA(rAU * self.UEqn.H(), self.U, self.p_rgh))

    #     # Buoyancy flux
    #     phig = surfaceScalarField(
    #         Word("phig"),
    #         -rAUf * self.ghf * fvc.snGrad(self.rhok) * self.mesh.magSf(),
    #     )

    #     # Flux based on HbyA, temporal correction, and buoyancy
    #     phiHbyA = surfaceScalarField(
    #         Word("phiHbyA"),
    #         fvc.flux(HbyA)
    #         + self.MRF.zeroFilter(rAUf * fvc.ddtCorr(self.U, self.phi))
    #         + phig,
    #     )

    #     self.MRF.makeRelative(phiHbyA)

    #     # Update pressure boundary conditions to ensure flux consistency
    #     constrainPressure(self.p_rgh, self.U, phiHbyA, rAUf, self.MRF)

    #     while pimple.correctNonOrthogonal():
    #         p_rghEqn = fvScalarMatrix(
    #             fvm.laplacian(rAUf, self.p_rgh) == fvc.div(phiHbyA)
    #         )

    #         p_rghEqn.setReference(self.pRefCell, self.pRefValue, False)
    #         p_rghEqn.solve(self.p_rgh.select(pimple.finalInnerIter()))

    #         if pimple.finalNonOrthogonalIter():
    #             # Conservative fluxes
    #             self.phi.assign(phiHbyA - p_rghEqn.flux())

    #     # Explicit relaxation for momentum corrector
    #     self.p_rgh.relax()

    #     # Correct momentum source with buoyancy & pressure correction
    #     self.U.assign(
    #         HbyA + rAU * fvc.reconstruct((phig - p_rghEqn.flux()) / rAUf)
    #     )
    #     self.U.correctBoundaryConditions()
    #     self.fvOptions.correct(self.U)

    #     # Correct face velocities if mesh is moving
    #     fvc.correctUf(self.Uf, self.U, self.phi)

    #     # Make fluxes relative to mesh motion
    #     fvc.makeRelative(self.phi, self.U)

    #     # Continuity errors (optional include)
    #     continuityErrs(self.phi, self.U, self.p_rgh)

    #     # Reconstruct physical pressure
    #     self.p.assign(self.p_rgh + self.rhok * self.gh)

    #     if self.p_rgh.needReference():
    #         self.p += dimensionedScalar(
    #             Word("p"),
    #             self.p.dimensions(),
    #             self.pRefValue - getRefCellValue(self.p, self.pRefCell),
    #         )
    #         self.p_rgh.assign(self.p - self.rhok * self.gh)
