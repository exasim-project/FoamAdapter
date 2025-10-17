from typing import Literal, Optional, Callable
from dataclasses import dataclass
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
from ...modules import fields

ControlDict = create_model(
    "controlDict",
    maxCo=(Optional[float], Field(default=None, description="Maximum Courant number")),
    adjustTimeStep=(
        Literal["yes", "no"],
        Field(default="no", description="Adjust time step"),
    ),
    __base__=ControlDictBase,
)

divSchemes = create_model(
    "divSchemes",
    __base__=DIVSchemes,
)

FvSchemes = create_model("fvSchemes", __base__=FvSchemesBase)

@dataclass
class PressureVelocityFields:
    p: volScalarField
    U: volVectorField
    phi: surfaceScalarField
    turbulence: any  # Placeholder for turbulence model

    # @staticmethod
    # def create(fields,models):
    #     p = volScalarField.read_field(fields.mesh, "p")
    #     U = volVectorField.read_field(fields.mesh, "U")
    #     phi = createPhi(U)
    #     turbulence = models.get("turbulence", None)
    #     return PressureVelocityFields(p=p, U=U, phi=phi, turbulence=turbulence)

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
    
    @staticmethod
    def register_fields(mesh: fvMesh) -> dict[str, Callable[[], fields.Field]]:
        required_fields = {
            "p": fields.FieldReader(field_class=volScalarField, mesh=mesh, field_name="p"),
            "U": fields.FieldReader(field_class=volVectorField, mesh=mesh, field_name="U"),
            "phi": fields.FieldReader(field_class=surfaceScalarField, mesh=mesh, field_name="phi"),
        }
        return required_fields

    def __init__(self, mesh, turbulence=None, laminarTransport=None):
        """
        Initialize the PIMPLE algorithm with the mesh and optional turbulence model.
        :param mesh: The computational mesh.
        :param turbulence: Optional turbulence model.
        :param laminarTransport: Optional laminar transport model.
        """
        self.mesh = mesh
        # self.turbulence = turbulence
        # self.laminarTransport = laminarTransport
        # self.p = volScalarField.read_field(mesh, "p")
        # self.U = volVectorField.read_field(mesh, "U")
        # self.phi = createPhi(self.U)
        self.UEqn = None  # Placeholder for momentum equation

        self.fvSolution = dictionary.read("system/fvSolution")
        self.pRefCell, self.pRefValue = None, None  # setRefCell(self.p, self.fvSolution.subDict("PIMPLE"))
        self.mesh.setFluxRequired(Word("p"))

    def stability_criteria(self, fields: PressureVelocityFields) -> CFLNumber:
        return CFLNumber(fields.phi, 5.0)

    def momentum_equation(self, pimple, fields: PressureVelocityFields) -> None:
        """
        Solve the momentum equations using the PIMPLE algorithm.
        """
        U, p, phi, turbulence = fields.U, fields.p, fields.phi, fields.turbulence
        if self.pRefCell is None or self.pRefValue is None:
            self.pRefCell, self.pRefValue = setRefCell(p, self.fvSolution.subDict("PIMPLE"))

        self.UEqn = fvVectorMatrix(
            fvm.ddt(U)
            + fvm.div(phi, U)
            + turbulence.divDevReff(U)
        )

        self.UEqn.relax()

        if pimple.momentumPredictor():
            solve(self.UEqn + fvc.grad(p))

    def pressure_correction(self, pimple, fields: PressureVelocityFields) -> None:
        """
        Correct the solution based on the PIMPLE algorithm.
        """
        U, p, phi = fields.U, fields.p, fields.phi
        rAU = volScalarField(Word("rAU"), 1.0 / self.UEqn.A())
        HbyA = volVectorField(constrainHbyA(rAU * self.UEqn.H(), U, p))

        phiHbyA = surfaceScalarField(
            Word("phiHbyA"),
            fvc.flux(HbyA) + fvc.interpolate(rAU) * fvc.ddtCorr(U, phi),
        )

        adjustPhi(phiHbyA, U, p)
        constrainPressure(p, U, phiHbyA, rAU)

        while pimple.correctNonOrthogonal():
            pEqn = fvScalarMatrix(fvm.laplacian(rAU, p) - fvc.div(phiHbyA))
            pEqn.setReference(self.pRefCell, self.pRefValue, False)
            pEqn.solve(p.select(pimple.finalInnerIter()))

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
