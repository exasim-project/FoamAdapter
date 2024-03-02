// SPDX-License-Identifier: GPLv-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM_GPL/conversion/convert.hpp"
#include "NeoFOAM/fields/field.hpp"

// namespace NeoFOAM
// {


template <typename NFoamType,typename FoamType>
NeoFOAM::Field<NFoamType> fromFoamField(const NeoFOAM::executor exec, const Foam::Field<FoamType> &field)
{
    // Create device-side views
    NeoFOAM::Field<NFoamType> nfField(exec, field.size());

    Kokkos::View<NFoamType *, Kokkos::HostSpace> CPU_view("copy_Host", field.size());
    for (Foam::label i = 0; i < field.size(); i++)
    {
        if constexpr (!std::is_same<FoamType, NFoamType>::value)
        {
            CPU_view(i) = convert(field[i]);
        }
        else
        {
            CPU_view(i) = field[i];
        }
    }
    if (std::holds_alternative<NeoFOAM::GPUExecutor>(exec))
    {
        Kokkos::View<NFoamType *, NeoFOAM::GPUExecutor::exec, Kokkos::MemoryUnmanaged>
            GPU_view(nfField.data(), field.size());
        Kokkos::deep_copy(GPU_view,CPU_view);
    }
    else
    {
        Kokkos::View<NFoamType *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
            N_CPU_view(nfField.data(), field.size());
        Kokkos::deep_copy(N_CPU_view,CPU_view);
    }

    return nfField;
};

template <typename NFoamType,typename FoamType>
NeoFOAM::Field<NFoamType> fromFoamField(const NeoFOAM::executor &exec, const Foam::List<FoamType> &field)
{
    // Create device-side views
    NeoFOAM::Field<NFoamType> nfField(exec, field.size());

    NeoFOAM::Field<NFoamType> CPUField(NeoFOAM::CPUExecutor {}, field.size());
    Kokkos::View<NFoamType *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> CPU_view(CPUField.data(), field.size());
    for (Foam::label i = 0; i < field.size(); i++)
    {
        if constexpr (!std::is_same<FoamType, NFoamType>::value)
        {
            CPU_view(i) = convert(field[i]);
        }
        else
        {
            CPU_view(i) = field[i];
        }
    }
    if (std::holds_alternative<NeoFOAM::GPUExecutor>(exec))
    {
        Kokkos::View<NFoamType *, NeoFOAM::GPUExecutor::exec, Kokkos::MemoryUnmanaged>
            GPU_view(nfField.data(), field.size());
        Kokkos::deep_copy(GPU_view,CPU_view);
    }
    else
    {
        Kokkos::View<NFoamType *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
            N_CPU_view(nfField.data(), field.size());
        Kokkos::deep_copy(N_CPU_view,CPU_view);
    }

    return nfField;
};
