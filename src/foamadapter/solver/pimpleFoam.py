import io
import pybFoam
import sys
import pybFoam
from pydantic import create_model, Field
from pybFoam.io.model_base import IOModelBase
from ..inputs_files.system import ControlDictBase, FvSchemesBase, DIVSchemes
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
    dynamicFvMesh,
    solve,
    adjustPhi,
    constrainPressure,
    createPhi,
    setRefCell,
    constrainHbyA,
    pimpleControl,
    computeCFLNumber,
)
from pybFoam.turbulence import incompressibleTurbulenceModel, singlePhaseTransportModel
from ..modules.pressureVelocityCoupling.incompressible import PimpleAlgorithm
from ..modules.stability_criteria import StabilityCriteria

ControlDict = create_model(
    "controlDict",
    maxCo=(float, Field(..., description="Maximum Courant number")),
    __base__=ControlDictBase,
)

divSchemes = create_model(
    "divSchemes",
    __base__=DIVSchemes,
)

FvSchemes = create_model("fvSchemes", __base__=FvSchemesBase)


class TransportProperties(IOModelBase):
    transportModel: str = Field(..., description="Transport model type")
    nu: float = Field(..., description="Kinematic viscosity", gt=0)


class _CaseInputs(IOModelBase):
    controlDict: ControlDict
    fvSchemes: FvSchemes
    transportProperties: TransportProperties


def create_fields(mesh):
    pU = PimpleAlgorithm(mesh=mesh)
    U = pU.U
    phi = pU.phi

    laminarTransport = singlePhaseTransportModel(U, phi)

    turbulence = incompressibleTurbulenceModel.New(U, phi, laminarTransport)
    pU.laminarTransport = laminarTransport
    pU.turbulence = turbulence

    return pU


class PimpleFoam:

    def inputs(self):
        """
        Return the input parameters for the PIMPLE algorithm.
        """
        return _CaseInputs

    def __init__(self, argv):
        """ """
        self.argList = pybFoam.argList(argv)

    def run(self):
        """
        Run the PIMPLE algorithm for the given mesh and turbulence model.
        :param pimple: The PIMPLE control object.
        """

        runTime = Time(self.argList)
        mesh = dynamicFvMesh.New(self.argList, runTime)

        pU = create_fields(mesh)
        stability_criteria = StabilityCriteria()
        stability_criteria.add_criteria(pU.stability_criteria())

        pimple = pimpleControl(mesh)

        while runTime.loop():
            Info(f"Time = {runTime.timeName()}")

            # Compute Courant number
            stability_criteria.setDelta(runTime)

            while pimple.loop():

                pU.momentum_equation(pimple)

                while pimple.correct():
                    pU.pressure_correction(pimple)

                if pimple.turbCorr():
                    pU.laminarTransport.correct()
                    pU.turbulence.correct()

            runTime.write(True)
            runTime.printExecutionTime()

        Info("End")


def main(argv):

    PimpleFoam(argv).run()


if __name__ == "__main__":
    main(sys.argv)
