// SPDX-License-Identifier: GPL-3.0-or-later
//
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
/* This file implements comparison operator to compare OpenFOAM and corresponding NeoFOAM fields
 * TODO the comparison operator only make sense for testing purposes
 * so this should be part of the tets
 */
#pragma once

#include <span>

#include "Field.H"

#include "NeoN/NeoN.hpp"

#include "FoamAdapter/conversion/convert.hpp"
#include "FoamAdapter/conversion/type_conversion.hpp"

namespace Foam
{

namespace fvcc = NeoN::finiteVolume::cellCentred;

#define FIELD_EQUALITY_OPERATOR(NF_TYPE, OF_TYPE)                                                  \
    bool operator==(const NeoN::Vector<NF_TYPE>& nf, const Foam::Field<OF_TYPE>& of)               \
    {                                                                                              \
        auto nfHost = nf.copyToHost();                                                             \
        auto nfView = nfHost.view();                                                               \
        for (int i = 0; i < nfView.size(); i++)                                                    \
        {                                                                                          \
            if (nfView[i] != convert(of[i]))                                                       \
            {                                                                                      \
                return false;                                                                      \
            }                                                                                      \
        }                                                                                          \
        return true;                                                                               \
    };


FIELD_EQUALITY_OPERATOR(NeoN::label, Foam::label)
FIELD_EQUALITY_OPERATOR(NeoN::scalar, Foam::scalar)
FIELD_EQUALITY_OPERATOR(NeoN::Vec3, Foam::vector)

#define VOLGEOFIELD_EQUALITY_OPERATOR(NF_TYPE, OF_TYPE)                                            \
    bool operator==(                                                                               \
        fvcc::VolumeField<NF_TYPE>& nf,                                                            \
        const Foam::GeometricField<OF_TYPE, Foam::fvPatchField, Foam::volMesh>& of                 \
    )                                                                                              \
    {                                                                                              \
        if (nf.internalVector() != of.internalField())                                             \
        {                                                                                          \
            return false;                                                                          \
        }                                                                                          \
                                                                                                   \
        /* compare boundaryVector */                                                               \
        /* NeoFOAM boundaries are stored in contiguous memory */                                   \
        /* whereas OpenFOAM boundaries are stored in a vector of patches */                        \
        auto nfBoundaryHost = nf.boundaryVector().value().copyToHost();                            \
        auto nfBoundaryView = nfBoundaryHost.view();                                               \
        NeoN::label pFacei = 0;                                                                    \
        for (const auto& patch : of.boundaryField())                                               \
        {                                                                                          \
            int patchSize = patch.size();                                                          \
            for (const auto& patchValue : patch)                                                   \
            {                                                                                      \
                if (nfBoundaryView[pFacei] != Foam::convert(patchValue))                           \
                {                                                                                  \
                    return false;                                                                  \
                }                                                                                  \
                pFacei++;                                                                          \
            }                                                                                      \
        }                                                                                          \
                                                                                                   \
        return true;                                                                               \
    }

VOLGEOFIELD_EQUALITY_OPERATOR(NeoN::scalar, Foam::scalar)
VOLGEOFIELD_EQUALITY_OPERATOR(NeoN::Vec3, Foam::vector)

/*
 *
 */
#define SURFGEOFIELD_EQUALITY_OPERATOR(NF_TYPE, OF_TYPE)                                           \
    bool operator==(                                                                               \
        const fvcc::SurfaceField<NF_TYPE>& nf,                                                     \
        const Foam::GeometricField<OF_TYPE, Foam::fvsPatchField, Foam::surfaceMesh>& of            \
    )                                                                                              \
    {                                                                                              \
        if (nf.internalVector() != of.internalField())                                             \
        {                                                                                          \
            return false;                                                                          \
        }                                                                                          \
                                                                                                   \
        /* compare boundaryVector */                                                               \
        /* NeoFOAM boundaries are stored in contiguous memory */                                   \
        /* whereas OpenFOAM boundaries are stored in a vector of patches */                        \
        auto nfHost = nf.internalVector().copyToHost();                                            \
        auto nfView = nfHost.view();                                                               \
        NeoN::label nInternalFaces = nf.internalVector().size();                                   \
        NeoN::label pFacei = nInternalFaces;                                                       \
        for (const auto& patch : of.boundaryField())                                               \
        {                                                                                          \
            int patchSize = patch.size();                                                          \
            for (const auto& patchValue : patch)                                                   \
            {                                                                                      \
                if (nfView[pFacei] != Foam::convert(patchValue))                                   \
                {                                                                                  \
                    return false;                                                                  \
                }                                                                                  \
                pFacei++;                                                                          \
            }                                                                                      \
        }                                                                                          \
                                                                                                   \
        return true;                                                                               \
    }

SURFGEOFIELD_EQUALITY_OPERATOR(NeoN::scalar, Foam::scalar)
SURFGEOFIELD_EQUALITY_OPERATOR(NeoN::Vec3, Foam::vector)

} // namespace Foam
