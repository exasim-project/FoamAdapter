// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
#pragma once

#include "volFields.H"
#include "surfaceFields.H"

#include "NeoN/NeoN.hpp"

#include "FoamAdapter/auxiliary/convert.hpp"

namespace FoamAdapter
{

namespace fvcc = NeoN::finiteVolume::cellCentred;
template<typename From>
struct type_map
{
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Foam::GeometricField<Foam::scalar, Foam::fvPatchField, Foam::volMesh>>
{
    using container_type = fvcc::VolumeField<NeoN::scalar>;
    using mapped_type = NeoN::scalar;
};

template<>
struct type_map<Foam::GeometricField< Foam::vector,  Foam::fvPatchField,  Foam::volMesh>>
{
    using container_type = fvcc::VolumeField<NeoN::Vec3>;
    using mapped_type = NeoN::Vec3;
};

template<>
struct type_map<Foam::GeometricField< Foam::scalar,  Foam::fvsPatchField,  Foam::surfaceMesh>>
{
    using container_type = fvcc::SurfaceField<NeoN::scalar>;
    using mapped_type = NeoN::scalar;
};

template<>
struct type_map<Foam::GeometricField< Foam::vector,  Foam::fvsPatchField,  Foam::surfaceMesh>>
{
    using container_type = fvcc::SurfaceField<NeoN::Vec3>;
    using mapped_type = NeoN::Vec3;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map< Foam::Field< Foam::scalar>>
{
    using container_type = NeoN::Vector<NeoN::scalar>;
    using mapped_type = NeoN::scalar;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map< Foam::Field<Foam::vector>>
{
    using container_type = NeoN::Vector<NeoN::Vec3>;
    using mapped_type = NeoN::Vec3;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Foam::List< Foam::scalar>>
{
    using container_type = NeoN::Vector<NeoN::scalar>;
    using mapped_type = NeoN::scalar;
};

// Specializations of type_map for specific type mappings.
template<>
struct type_map<Foam::List<Foam::label>>
{
    using container_type = NeoN::Vector<NeoN::label>;
    using mapped_type = NeoN::label;
};

}; // namespace Foam
