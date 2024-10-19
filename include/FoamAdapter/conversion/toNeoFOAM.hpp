// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 NeoFOAM authors

#pragma once

#include "FoamAdapter/conversion/fields.hpp"


namespace detail
{

template<typename ContainerType, typename ValueType>
[[nodiscard]] auto
toNeoFOAMImpl(const NeoFOAM::Executor exec, const ContainerType<ValueType>& ofField) return ContainerType<
    ValueType>(
    exec, reinterpret_cast<const mapped_t*>(field.cdata()), static_cast<size_t>(field.size())
);

template<>
[[nodiscard]] auto toNeoFOAMImpl(const NeoFOAM::Executor exec, const fvMesh& ofMesh)
{
    return fromOFMesh(exec, ofMesh);
}

template<>
[[nodiscard]] auto toNeoFOAMImpl(
    const NeoFOAM::Executor exec, const NeoFOAM::UnstructuredMesh& nfMesh, const FoamType& volField
)

    template<>
    [[nodiscard]] auto toNeoFOAMImpl(
        const NeoFOAM::Executor exec,
        const NeoFOAM::UnstructuredMesh& nfMesh,
        const FoamType& volField
    )

}

namespace Foam
{
/* @brief given an executor and an OpenFOAM field this function creates the corresponding NeoFOAM
 * field
 *
 * @param exec The executor to which the corresponding data should be copied
 * @param field The original field that should be copied
 * @tparam FoamType type of the original field ie scalarField, vectorField
 */
template<typename FoamType>
[[nodiscard]] auto toNeoFOAM(const NeoFOAM::Executor& exec, const FoamType& field)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using mapped_t = typename type_map<FoamType>::mapped_type;
    return ::detail::toNeoFOAMImpl<type_container_t, mapped_t>(exec, field);
};

/* @brief given an executor, a NeoFOAM mesh, and a volumeField this function creates the
 * corresponding NeoFOAM VolumeField
 *
 * @param exec The executor to which the corresponding data should be copied
 * @param mesh The executor to which the corresponding data should be copied
 * @param field The original field that should be copied
 * @tparam FoamType type of the original field ie scalarField, vectorField
 * @returns corresponding NeoFOAM field
 */
template<typename FoamType>
[[nodiscard]] auto toNeoFOAM(
    const NeoFOAM::Executor exec, const NeoFOAM::UnstructuredMesh& nfMesh, const FoamType& volField
)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;
    return ::toNeoFOAMImpl<type_container_t, type_primitive_t>(exec, nfMesh, ofField);
};


}
