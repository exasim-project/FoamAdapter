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

namespace fvcc = NeoN::finiteVolume::cellCentred;
template<typename From>
struct type_map
{
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<GeometricField<scalar, fvPatchField, volMesh>>
{
    using container_type = fvcc::VolumeField<NeoN::scalar>;
    using mapped_type = NeoN::scalar;
};

template<>
struct type_map<GeometricField<vector, fvPatchField, volMesh>>
{
    using container_type = fvcc::VolumeField<NeoN::Vec3>;
    using mapped_type = NeoN::Vec3;
};

template<>
struct type_map<GeometricField<scalar, fvsPatchField, surfaceMesh>>
{
    using container_type = fvcc::SurfaceField<NeoN::scalar>;
    using mapped_type = NeoN::scalar;
};

template<>
struct type_map<GeometricField<vector, fvsPatchField, surfaceMesh>>
{
    using container_type = fvcc::SurfaceField<NeoN::Vec3>;
    using mapped_type = NeoN::Vec3;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Field<scalar>>
{
    using container_type = NeoN::Field<NeoN::scalar>;
    using mapped_type = NeoN::scalar;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Field<Vec3>>
{
    using container_type = NeoN::Field<NeoN::Vec3>;
    using mapped_type = NeoN::Vec3;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<List<scalar>>
{
    using container_type = NeoN::Field<NeoN::scalar>;
    using mapped_type = NeoN::scalar;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<List<label>>
{
    using container_type = NeoN::Field<NeoN::label>;
    using mapped_type = NeoN::label;
};

}; // namespace Foam
