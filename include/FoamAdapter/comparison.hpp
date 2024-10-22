// SPDX-License-Identifier: GPL-3.0-or-later
//
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
/* This file implements comparison operator to compare OpenFOAM and corresponding NeoFOAM fields
 */
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

#define FIELD_EQUALITY_OPERATOR(NF_TYPE, OF_TYPE)                                                  \
    bool operator==(const NeoFOAM::Field<NF_TYPE>& nf, const Foam::Field<OF_TYPE>& of)             \
    {                                                                                              \
        auto nfHost = nf.copyToHost();                                                             \
        auto nfSpan = nfHost.span();                                                               \
        for (int i = 0; i < nfSpan.size(); i++)                                                    \
        {                                                                                          \
            if (nfSpan[i] != Foam::convert(of[i]))                                                 \
            {                                                                                      \
                return false;                                                                      \
            }                                                                                      \
        }                                                                                          \
        return true;                                                                               \
    };


FIELD_EQUALITY_OPERATOR(NeoFOAM::scalar, Foam::scalar)
FIELD_EQUALITY_OPERATOR(NeoFOAM::Vector, Foam::vector)

#define VOLGEOFIELD_EQUALITY_OPERATOR(NF_TYPE, OF_TYPE)                                            \
    bool operator==(                                                                               \
        fvcc::VolumeField<NF_TYPE>& nf,                                                            \
        const Foam::GeometricField<OF_TYPE, Foam::fvPatchField, Foam::volMesh>& of                 \
    )                                                                                              \
    {                                                                                              \
        if (nf.internalField() != of.internalField())                                              \
        {                                                                                          \
            return false;                                                                          \
        }                                                                                          \
                                                                                                   \
        /* compare boundaryField */                                                                \
        /* NeoFOAM boundaries are stored in contiguous memory */                                   \
        /* whereas OpenFOAM boundaries are stored in a vector of patches */                        \
        auto nfBoundaryHost = nf.boundaryField().value().copyToHost();                             \
        auto nfBoundarySpan = nfBoundaryHost.span();                                               \
        NeoFOAM::label pFacei = 0;                                                                 \
        for (const auto& patch : of.boundaryField())                                               \
        {                                                                                          \
            int patchSize = patch.size();                                                          \
            for (const auto& patchValue : patch)                                                   \
            {                                                                                      \
                if (nfBoundarySpan[pFacei] != Foam::convert(patchValue))                           \
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

/*
 *
 */
#define SURFGEOFIELD_EQUALITY_OPERATOR(NF_TYPE, OF_TYPE)                                           \
    bool operator==(                                                                               \
        const fvcc::SurfaceField<NF_TYPE>& nf,                                                     \
        const Foam::GeometricField<OF_TYPE, Foam::fvsPatchField, Foam::surfaceMesh>& of            \
    )                                                                                              \
    {                                                                                              \
        if (nf.internalField() != of.internalField())                                              \
        {                                                                                          \
            return false;                                                                          \
        }                                                                                          \
                                                                                                   \
        /* compare boundaryField */                                                                \
        /* NeoFOAM boundaries are stored in contiguous memory */                                   \
        /* whereas OpenFOAM boundaries are stored in a vector of patches */                        \
        auto nfHost = nf.internalField().copyToHost();                                             \
        auto nfSpan = nfHost.span();                                                               \
        NeoFOAM::label nInternalFaces = nf.internalField().size();                                 \
        NeoFOAM::label pFacei = nInternalFaces;                                                    \
        for (const auto& patch : of.boundaryField())                                               \
        {                                                                                          \
            int patchSize = patch.size();                                                          \
            for (const auto& patchValue : patch)                                                   \
            {                                                                                      \
                if (nfSpan[pFacei] != Foam::convert(patchValue))                                   \
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
