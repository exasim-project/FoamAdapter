// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/fvcc/surfaceInterpolation/surfaceInterpolationFactory.hpp"

namespace Foam
{

// defineTypeNameAndDebug(surfaceInterpolationFactory, 0);
defineRunTimeSelectionTable(surfaceInterpolationFactory, dictionary);

} // namespace NeoFOAM

Foam::autoPtr<Foam::SurfaceInterpolationKernel> Foam::surfaceInterpolationFactory::New(
    const NeoFOAM::Executor& exec, const NeoFOAM::UnstructuredMesh& mesh
)
{
    word schemeType("isoAlpha");

    Info << "Selecting SurfaceInterpolationKernel: " << schemeType << endl;

    auto* ctorPtr = dictionaryConstructorTable(schemeType);

    if (!ctorPtr)
    {
        FatalErrorInLookup(
            "SurfaceInterpolationKernels", schemeType, *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<SurfaceInterpolationKernel>(ctorPtr(exec, mesh));
}
