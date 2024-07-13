// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "FoamAdapter/conversion/convert.hpp"
#include "NeoFOAM/finiteVolume/cellCentred.hpp"
#include "NeoFOAM/core/dictionary.hpp"
#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/field.hpp"

#include "FoamAdapter/conversion/type_conversion.hpp"

namespace Foam
{
namespace fvcc = NeoFOAM::finiteVolume::cellCentred;
template<typename FoamType>
auto fromFoamField(const NeoFOAM::Executor& exec, const FoamType& field)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;
    // Create device-side views
    type_container_t nfField(exec, field.size());

    type_container_t CPUField(NeoFOAM::CPUExecutor {}, field.size());
    Kokkos::View<type_primitive_t*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> CPU_view(
        CPUField.data(), field.size()
    );
    for (Foam::label i = 0; i < field.size(); i++)
    {
        CPU_view(i) = convert(field[i]);
    }
    if (std::holds_alternative<NeoFOAM::GPUExecutor>(exec))
    {
        Kokkos::View<type_primitive_t*, NeoFOAM::GPUExecutor::exec, Kokkos::MemoryUnmanaged>
            GPU_view(nfField.data(), field.size());
        Kokkos::deep_copy(GPU_view, CPU_view);
    }
    else
    {
        Kokkos::View<type_primitive_t*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> N_CPU_view(
            nfField.data(), field.size()
        );
        Kokkos::deep_copy(N_CPU_view, CPU_view);
    }

    return nfField;
};

template<typename FoamType>
auto readVolBoundaryConditions(const NeoFOAM::UnstructuredMesh& uMesh, const FoamType& volField)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;

    std::vector<std::unique_ptr<fvcc::VolumeBoundary<type_primitive_t>>> bcs;

    // get boundary as dictionary
    Foam::OStringStream os;
    volField.boundaryField().writeEntries(os);
    Foam::IStringStream is(os.str());
    Foam::dictionary bDict(is);
    Foam::wordList bNames = bDict.toc();
    int patchi = 0;

    for (const auto& bName : bNames)
    {
        Foam::Info << "Boundary name: " << bName << Foam::endl;
        Foam::dictionary patchDict = bDict.subDict(bName);
        Foam::Info << "Boundary type: " << patchDict.get<Foam::word>("type") << Foam::endl;
        Foam::word type = patchDict.get<Foam::word>("type");
        NeoFOAM::Dictionary neoPatchDict;
        neoPatchDict.insert("type", std::string(type));
        bcs.push_back(VolumeBoundary<type_primitive_t>(uMesh, patchi, neoPatchDict));
        patchi++;
    }
    return bcs;
}

template<typename FoamType>
auto constructFrom(
    const NeoFOAM::Executor exec, const NeoFOAM::UnstructuredMesh& uMesh, const FoamType& volField
)
{

    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;

    type_container_t nfVolField(exec, uMesh, std::move(readVolBoundaryConditions(uMesh, volField)));

    nfVolField.internalField() = fromFoamField(exec, volField.primitiveField());
    nfVolField.correctBoundaryConditions();

    return nfVolField;
};

template<typename FoamType>
auto readSurfaceBoundaryConditions(
    const NeoFOAM::UnstructuredMesh& uMesh, const FoamType& surfaceField
)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;

    std::vector<std::unique_ptr<fvcc::SurfaceBoundary<type_primitive_t>>> bcs;

    // get boundary as dictionary
    Foam::OStringStream os;
    surfaceField.boundaryField().writeEntries(os);
    Foam::IStringStream is(os.str());
    Foam::dictionary bDict(is);
    Foam::wordList bNames = bDict.toc();
    int patchi = 0;

    for (const auto& bName : bNames)
    {
        Foam::Info << "Boundary name: " << bName << Foam::endl;
        Foam::dictionary patchDict = bDict.subDict(bName);
        Foam::Info << "Boundary type: " << patchDict.get<Foam::word>("type") << Foam::endl;
        Foam::word type = patchDict.get<Foam::word>("type");
        NeoFOAM::Dictionary npatchDict;
        npatchDict.insert("type", std::string(type));
        bcs.push_back(SurfaceBoundary<type_primitive_t>(uMesh, patchi, npatchDict));
        patchi++;
    }
    return bcs;
}

template<typename FoamType>
auto constructSurfaceField(
    const NeoFOAM::Executor exec, const NeoFOAM::UnstructuredMesh& uMesh, const FoamType& surfField
)
{

    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;
    using foam_primitive_t = typename FoamType::cmptType;

    type_container_t nfSurfField(
        exec, uMesh, std::move(readSurfaceBoundaryConditions(uMesh, surfField))
    );

    Field<foam_primitive_t> flattenedField(nfSurfField.internalField().size());
    size_t nInternal = uMesh.nInternalFaces();

    forAll(surfField, facei)
    {
        flattenedField[facei] = convert(surfField[facei]);
    }

    Foam::label idx = nInternal;
    forAll(surfField.boundaryField(), patchi)
    {
        const fvsPatchField<foam_primitive_t>& psurfField = surfField.boundaryField()[patchi];

        forAll(psurfField, facei)
        {
            flattenedField[idx] = psurfField[facei];
            idx++;
        }
    }
    assert(idx == flattenedField.size());

    nfSurfField.internalField() = fromFoamField(exec, flattenedField);
    nfSurfField.correctBoundaryConditions();

    return nfSurfField;
}

}; // namespace Foam
