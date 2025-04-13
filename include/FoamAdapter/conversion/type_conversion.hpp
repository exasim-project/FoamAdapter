// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "volFields.H"
#include "surfaceFields.H"

#include "NeoN/finiteVolume/cellCentred/fields/volumeField.hpp"
#include "NeoN/finiteVolume/cellCentred/fields/surfaceField.hpp"
#include "NeoN/core/executor/executor.hpp"
#include "NeoN/fields/field.hpp"

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
struct type_map<GeometricField<scalar, fvPatchField, volMesh>>
{
    using container_type = fvcc::VolumeField<NeoFOAM::scalar>;
    using mapped_type = NeoFOAM::scalar;
};

template<>
struct type_map<GeometricField<vector, fvPatchField, volMesh>>
{
    using container_type = fvcc::VolumeField<NeoFOAM::Vector>;
    using mapped_type = NeoFOAM::Vector;
};

template<>
struct type_map<GeometricField<scalar, fvsPatchField, surfaceMesh>>
{
    using container_type = fvcc::SurfaceField<NeoFOAM::scalar>;
    using mapped_type = NeoFOAM::scalar;
};

template<>
struct type_map<GeometricField<vector, fvsPatchField, surfaceMesh>>
{
    using container_type = fvcc::SurfaceField<NeoFOAM::Vector>;
    using mapped_type = NeoFOAM::Vector;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Field<scalar>>
{
    using container_type = NeoFOAM::Field<NeoFOAM::scalar>;
    using mapped_type = NeoFOAM::scalar;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Field<vector>>
{
    using container_type = NeoFOAM::Field<NeoFOAM::Vector>;
    using mapped_type = NeoFOAM::Vector;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<List<scalar>>
{
    using container_type = NeoFOAM::Field<NeoFOAM::scalar>;
    using mapped_type = NeoFOAM::scalar;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<List<label>>
{
    using container_type = NeoFOAM::Field<NeoFOAM::label>;
    using mapped_type = NeoFOAM::label;
};

}; // namespace Foam
