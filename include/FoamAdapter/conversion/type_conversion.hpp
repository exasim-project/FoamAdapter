// SPDX-License-Identifier: GPLv-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "FoamAdapter/conversion/convert.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccVolField.hpp"
#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"
#include "volFields.H"

namespace Foam {
template <typename From> struct type_map {};

// Specializations of type_map for specific type mappings.
template <>
struct type_map<
    Foam::GeometricField<Foam::scalar, Foam::fvPatchField, Foam::volMesh>> {
  using container_type = NeoFOAM::fvccVolField<NeoFOAM::scalar>;
  using mapped_type = NeoFOAM::scalar;
};

template <>
struct type_map<
    Foam::GeometricField<Foam::vector, Foam::fvPatchField, Foam::volMesh>> {
  using container_type = NeoFOAM::fvccVolField<NeoFOAM::Vector>;
  using mapped_type = NeoFOAM::Vector;
};

// Specializations of type_map for specific type mappings.
template <> struct type_map<Foam::Field<Foam::scalar>> {
  using container_type = NeoFOAM::Field<NeoFOAM::scalar>;
  using mapped_type = NeoFOAM::scalar;
};

// Specializations of type_map for specific type mappings.
template <> struct type_map<Foam::Field<Foam::vector>> {
  using container_type = NeoFOAM::Field<NeoFOAM::Vector>;
  using mapped_type = NeoFOAM::Vector;
};

// Specializations of type_map for specific type mappings.
template <> struct type_map<Foam::List<Foam::scalar>> {
  using container_type = NeoFOAM::Field<NeoFOAM::scalar>;
  using mapped_type = NeoFOAM::scalar;
};

// Specializations of type_map for specific type mappings.
template <> struct type_map<Foam::List<Foam::label>> {
  using container_type = NeoFOAM::Field<NeoFOAM::label>;
  using mapped_type = NeoFOAM::label;
};

}; // namespace Foam