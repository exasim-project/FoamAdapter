// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "volFields.H"
#include "surfaceFields.H"

#include "NeoFOAM/finiteVolume/cellCentred.hpp"
#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/field.hpp"

#include "FoamAdapter/conversion/convert.hpp"

namespace Foam
{

namespace fvcc = NeoFOAM::finiteVolume::cellCentred;
template<typename From>
struct type_map
{
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Foam::GeometricField<Foam::scalar, Foam::fvPatchField, Foam::volMesh>>
{
    using container_type = fvcc::VolumeField<NeoFOAM::scalar>;
    using mapped_type = NeoFOAM::scalar;
};

template<>
struct type_map<Foam::GeometricField<Foam::vector, Foam::fvPatchField, Foam::volMesh>>
{
    using container_type = fvcc::VolumeField<NeoFOAM::Vector>;
    using mapped_type = NeoFOAM::Vector;
};

template<>
struct type_map<Foam::GeometricField<Foam::scalar, Foam::fvsPatchField, Foam::surfaceMesh>>
{
    using container_type = fvcc::SurfaceField<NeoFOAM::scalar>;
    using mapped_type = NeoFOAM::scalar;
};

template<>
struct type_map<Foam::GeometricField<Foam::vector, Foam::fvsPatchField, Foam::surfaceMesh>>
{
    using container_type = fvcc::SurfaceField<NeoFOAM::Vector>;
    using mapped_type = NeoFOAM::Vector;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Foam::Field<Foam::scalar>>
{
    using container_type = NeoFOAM::Field<NeoFOAM::scalar>;
    using mapped_type = NeoFOAM::scalar;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Foam::Field<Foam::vector>>
{
    using container_type = NeoFOAM::Field<NeoFOAM::Vector>;
    using mapped_type = NeoFOAM::Vector;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Foam::List<Foam::scalar>>
{
    using container_type = NeoFOAM::Field<NeoFOAM::scalar>;
    using mapped_type = NeoFOAM::scalar;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Foam::List<Foam::label>>
{
    using container_type = NeoFOAM::Field<NeoFOAM::label>;
    using mapped_type = NeoFOAM::label;
};

}; // namespace Foam
