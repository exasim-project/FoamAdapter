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
from ..inputs_files.case_inputs import Registry, FileSpec
from ..turbulence.incompressible import TurbulenceModel
from ..modules.pressureVelocityCoupling.incompressible import PimpleAlgorithm
from ..modules.stability_criteria import StabilityCriteria
from ..modules.transportModels import singlePhaseTransportModel

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




def create_fields(mesh):
    pU = PimpleAlgorithm(mesh=mesh)
    U = pU.U
    phi = pU.phi

    transportModel = singlePhaseTransportModel(U, phi)

    turbulence = TurbulenceModel.New(U, phi, transportModel)

    pU.laminarTransport = transportModel
    pU.turbulence = turbulence

    return pU


class PimpleFoam:

    def inputs(self):
        """
        Return the input parameters for the PIMPLE algorithm.
        """
        registry = Registry({ })
        registry = PimpleAlgorithm.inputs(registry)
        registry = TurbulenceModel.inputs(registry)
        registry = singlePhaseTransportModel.inputs(registry)
        return registry

    def __init__(self, argv):
        """ """
        self._argv = argv

    def run(self):
        """
        Run the PIMPLE algorithm for the given mesh and turbulence model.
        :param pimple: The PIMPLE control object.
        """

        argList = pybFoam.argList(self._argv)
        runTime = Time(argList)
        mesh = dynamicFvMesh.New(argList, runTime)

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
