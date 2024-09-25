// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include <span>

#include "Field.H"

#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/core/primitives/label.hpp"

#include "FoamAdapter/conversion/convert.hpp"
#include "FoamAdapter/conversion/type_conversion.hpp"

namespace Foam
{

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;

#define FIELD_EQUALITY_OPERATOR(NFIELD_TYPE, FFIELD_TYPE)                                          \
    bool operator==(NeoFOAM::Field<NFIELD_TYPE>& nfield, const Foam::Field<FFIELD_TYPE>& ffield)   \
    {                                                                                              \
        auto hostSpan = nfield.copyToHost();                                                       \
        auto span = hostSpan.span();                                                               \
                                                                                                   \
        for (int i = 0; i < span.size(); i++)                                                      \
        {                                                                                          \
            if (span[i] != Foam::convert(ffield[i]))                                               \
            {                                                                                      \
                return false;                                                                      \
            }                                                                                      \
        }                                                                                          \
        return true;                                                                               \
    };

FIELD_EQUALITY_OPERATOR(NeoFOAM::scalar, Foam::scalar)
FIELD_EQUALITY_OPERATOR(NeoFOAM::Vector, Foam::vector)

#define VOLGEOFIELD_EQUALITY_OPERATOR(NFIELD_TYPE, FFIELD_TYPE)                                    \
    bool operator==(                                                                               \
        fvcc::VolumeField<NFIELD_TYPE>& nfield,                                                    \
        const Foam::GeometricField<FFIELD_TYPE, Foam::fvPatchField, Foam::volMesh>& ffield         \
    )                                                                                              \
    {                                                                                              \
        auto hostSpan = nfield.internalField().copyToHost();                                       \
        auto span = hostSpan.span();                                                               \
        const auto& internalFField = ffield.internalField();                                       \
        /* compare internalField*/                                                                 \
        for (int i = 0; i < span.size(); i++)                                                      \
        {                                                                                          \
            if (span[i] != Foam::convert(internalFField[i]))                                       \
            {                                                                                      \
                return false;                                                                      \
            }                                                                                      \
        }                                                                                          \
        /* compare boundaryField */                                                                \
        /* NeoFOAM boundaries are stored in contiguous memory */                                   \
        /* whereas OpenFOAM boundaries are stored in a vector of patches */                        \
        auto patchValueHost = nfield.boundaryField().value().copyToHost();                         \
        auto patchValueSpan = patchValueHost.span();                                               \
        NeoFOAM::label pFacei = 0;                                                                 \
        for (const auto& patch : ffield.boundaryField())                                           \
        {                                                                                          \
            int patchSize = patch.size();                                                          \
            for (const auto& patchValue : patch)                                                   \
            {                                                                                      \
                if (patchValueSpan[pFacei] != Foam::convert(patchValue))                           \
                {                                                                                  \
                    return false;                                                                  \
                }                                                                                  \
                pFacei++;                                                                          \
            }                                                                                      \
        }                                                                                          \
                                                                                                   \
        return true;                                                                               \
    }

VOLGEOFIELD_EQUALITY_OPERATOR(NeoFOAM::scalar, Foam::scalar)
VOLGEOFIELD_EQUALITY_OPERATOR(NeoFOAM::Vector, Foam::vector)

#define SURFGEOFIELD_EQUALITY_OPERATOR(NFIELD_TYPE, FFIELD_TYPE)                                   \
    bool operator==(                                                                               \
        fvcc::SurfaceField<NFIELD_TYPE>& nfield,                                                   \
        const Foam::GeometricField<FFIELD_TYPE, Foam::fvsPatchField, Foam::surfaceMesh>& ffield    \
    )                                                                                              \
    {                                                                                              \
        auto hostSpan = nfield.internalField().copyToHost();                                       \
        auto span = hostSpan.span();                                                               \
        const auto& internalFField = ffield.internalField();                                       \
        /* compare internalField the fvccSurfaceField contains the boundaryValues                  \
         */                                                                                        \
        for (int i = 0; i < internalFField.size(); i++)                                            \
        {                                                                                          \
            if (span[i] != Foam::convert(internalFField[i]))                                       \
            {                                                                                      \
                return false;                                                                      \
            }                                                                                      \
        }                                                                                          \
        /* compare boundaryField */                                                                \
        /* NeoFOAM boundaries are stored in contiguous memory */                                   \
        /* whereas OpenFOAM boundaries are stored in a vector of patches */                        \
        NeoFOAM::label nInternalFaces = nfield.internalField().size();                             \
        NeoFOAM::label pFacei = nInternalFaces;                                                    \
        for (const auto& patch : ffield.boundaryField())                                           \
        {                                                                                          \
            int patchSize = patch.size();                                                          \
            for (const auto& patchValue : patch)                                                   \
            {                                                                                      \
                if (hostSpan[pFacei] != Foam::convert(patchValue))                                 \
                {                                                                                  \
                    return false;                                                                  \
                }                                                                                  \
                pFacei++;                                                                          \
            }                                                                                      \
        }                                                                                          \
                                                                                                   \
        return true;                                                                               \
    }

SURFGEOFIELD_EQUALITY_OPERATOR(NeoFOAM::scalar, Foam::scalar)
SURFGEOFIELD_EQUALITY_OPERATOR(NeoFOAM::Vector, Foam::vector)

} // namespace Foam
