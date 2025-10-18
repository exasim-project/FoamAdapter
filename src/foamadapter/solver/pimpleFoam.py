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
from ..modules.pressureVelocityCoupling.incompressible import (
    PimpleAlgorithm,
    PressureVelocityFields,
)
from ..modules.stability_criteria import StabilityCriteria
from ..modules.transportModels import SinglePhaseTransportModel
from ..modules.fields import Fields
from ..modules.models import Models
from ..modules.setup import initialize_containers
from ..modules.setup import visualize_dag, containers_to_deps

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


def create_fields_models(mesh):
    pU = PimpleAlgorithm(mesh=mesh)

    fields = Fields()
    fields.add_fields(pU.register_fields(mesh=mesh))

    models = Models(models={})
    models.add_model("pU", pU)

    singlePhaseTransportModel = SinglePhaseTransportModel()
    models.add_model("singlePhaseTransportModel", singlePhaseTransportModel)

    turbulence = TurbulenceModel()
    models.add_model("turbulence", turbulence)

    # visualize_dag(
    #     containers_to_deps(fields, models), "pimpleFoam_fields_models_dag", show=True
    # )
    initialize_containers(fields, models)

    return fields, models


class PimpleFoam:

    def inputs(self):
        """
        Return the input parameters for the PIMPLE algorithm.
        """
        registry = Registry({})
        registry = PimpleAlgorithm.inputs(registry)
        registry = TurbulenceModel.inputs(registry)
        registry = SinglePhaseTransportModel.inputs(registry)
        return registry

    def __init__(self, argv):
        """
        Initialize the PIMPLE solver with command line arguments.
        """
        self._argv = argv

    def run(self):
        """
        Run the PIMPLE algorithm for the given mesh and turbulence model.
        :param pimple: The PIMPLE control object.
        """

        argList = pybFoam.argList(self._argv)
        runTime = Time(argList)
        mesh = dynamicFvMesh.New(argList, runTime)

        fields, models = create_fields_models(mesh)
        print("Fields and models created.")
        pU: PimpleAlgorithm = models["pU"]

        puData = PressureVelocityFields(
            p=fields["p"],
            U=fields["U"],
            phi=fields["phi"],
            turbulence=models["turbulence"],
        )

        stability_criteria = StabilityCriteria()
        stability_criteria.add_criteria(pU.stability_criteria(puData))

        pimple = pimpleControl(mesh)

        while runTime.loop():
            Info(f"Time = {runTime.timeName()}")

            # Compute Courant number
            stability_criteria.setDelta(runTime)

            while pimple.loop():

                pU.momentum_equation(pimple, puData)

                while pimple.correct():
                    pU.pressure_correction(pimple, puData)

                if pimple.turbCorr():
                    models["singlePhaseTransportModel"].correct()
                    models["turbulence"].correct()

            runTime.write(True)
            runTime.printExecutionTime()

        Info("End")


def main(argv):

    PimpleFoam(argv).run()


if __name__ == "__main__":
    main(sys.argv)
