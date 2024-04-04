// SPDX-License-Identifier: MPL-2.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/fvcc/surfaceInterpolation/surfaceInterpolationFactory.hpp"

namespace Foam
{

// defineTypeNameAndDebug(surfaceInterpolationFactory, 0);
defineRunTimeSelectionTable(surfaceInterpolationFactory, dictionary);

} // namespace NeoFOAM

Foam::autoPtr<Foam::surfaceInterpolationKernel>
Foam::surfaceInterpolationFactory::New
(
    const NeoFOAM::executor& exec,
    const NeoFOAM::unstructuredMesh& mesh
)
{
    word schemeType("isoAlpha");

    Info<< "Selecting surfaceInterpolationKernel: " << schemeType << endl;

    auto* ctorPtr = dictionaryConstructorTable(schemeType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "surfaceInterpolationKernels",
            schemeType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<surfaceInterpolationKernel>(ctorPtr(exec,mesh));
}