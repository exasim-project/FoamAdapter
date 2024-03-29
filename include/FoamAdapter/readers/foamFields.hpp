// SPDX-License-Identifier: GPLv-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/core/executor/executor.hpp"
#include "FoamAdapter/conversion/convert.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/fields/fvccVolField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/fvccBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/scalar/fvccScalarFixedValueBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/scalar/fvccScalarZeroGradientBoundaryField.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/bcFields/scalar/fvccScalarEmptyBoundaryField.hpp"

#include "FoamAdapter/conversion/type_conversion.hpp"
// namespace NeoFOAM
// {

template <typename FoamType>
auto fromFoamField(const NeoFOAM::executor &exec, const FoamType& field)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitve_t = typename type_map<FoamType>::mapped_type;
    // Create device-side views
    type_container_t nfField(exec, field.size());

    type_container_t CPUField(NeoFOAM::CPUExecutor{}, field.size());
    Kokkos::View<type_primitve_t *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> CPU_view(CPUField.data(), field.size());
    for (Foam::label i = 0; i < field.size(); i++)
    {
        CPU_view(i) = convert(field[i]);
    }
    if (std::holds_alternative<NeoFOAM::GPUExecutor>(exec))
    {
        Kokkos::View<type_primitve_t *, NeoFOAM::GPUExecutor::exec, Kokkos::MemoryUnmanaged>
            GPU_view(nfField.data(), field.size());
        Kokkos::deep_copy(GPU_view, CPU_view);
    }
    else
    {
        Kokkos::View<type_primitve_t *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
            N_CPU_view(nfField.data(), field.size());
        Kokkos::deep_copy(N_CPU_view, CPU_view);
    }

    return nfField;
};

template <typename FoamType>
auto readBoundaryCondition(
    const NeoFOAM::unstructuredMesh& uMesh,
    const FoamType& volField)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitve_t = typename type_map<FoamType>::mapped_type;

    std::vector<std::unique_ptr<NeoFOAM::fvccBoundaryField<type_primitve_t>>> bcs;

    // get boundary as dictionary
    Foam::OStringStream os;
    volField.boundaryField().writeEntries(os);
    Foam::IStringStream is(os.str());
    Foam::dictionary bDict(is);
    Foam::wordList bNames = bDict.toc();
    int patchi = 0;

    for (const auto &bName : bNames)
    {
        Foam::Info << "Boundary name: " << bName << Foam::endl;
        Foam::dictionary patchDict = bDict.subDict(bName);
        Foam::Info << "Boundary type: " << patchDict.get<Foam::word>("type") << Foam::endl;
        Foam::word type = patchDict.get<Foam::word>("type");
        if (type == "zeroGradient")
        {
            bcs.push_back(std::make_unique<NeoFOAM::fvccScalarZeroGradientBoundaryField>(uMesh, patchi));
        }
        else if (type == "fixedValue")
        {
            bcs.push_back(std::make_unique<NeoFOAM::fvccScalarFixedValueBoundaryField>(uMesh, patchi, patchDict.get<Foam::scalar>("value")));
        }
        else if (type == "empty")
        {
            bcs.push_back(std::make_unique<NeoFOAM::fvccScalarEmptyBoundaryField>(uMesh, patchi));
        }
        else
        {
            Foam::FatalError << "Boundary type not supported" << Foam::endl;
        }
        patchi++;
    }
    return bcs;
}


template <typename FoamType>
auto constructFrom(const NeoFOAM::executor exec, const NeoFOAM::unstructuredMesh &uMesh, const FoamType &volField)
{

    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitve_t = typename type_map<FoamType>::mapped_type;


    type_container_t nfVolField(
        exec,
        uMesh,
        std::move(readBoundaryCondition(uMesh,volField))
    );

    nfVolField.internalField() = fromFoamField(exec, volField.primitiveField());

    return nfVolField;
};