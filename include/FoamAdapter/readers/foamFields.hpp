// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/finiteVolume/cellCentred.hpp"
#include "NeoFOAM/core/dictionary.hpp"
#include "NeoFOAM/core/executor/executor.hpp"
#include "NeoFOAM/fields/field.hpp"

#include "FoamAdapter/conversion/convert.hpp"
#include "FoamAdapter/conversion/type_conversion.hpp"

namespace Foam
{
namespace fvcc = NeoFOAM::finiteVolume::cellCentred;
template<typename FoamType>
auto fromFoamField(const NeoFOAM::Executor& exec, const FoamType& field)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using mapped_t = typename type_map<FoamType>::mapped_type;
    type_container_t nfField(
        exec, reinterpret_cast<const mapped_t*>(field.cdata()), static_cast<size_t>(field.size())
    );

    return nfField;
};

template<typename FoamType>
auto readVolBoundaryConditions(const NeoFOAM::UnstructuredMesh& uMesh, const FoamType& volField)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;

    std::vector<fvcc::VolumeBoundary<type_primitive_t>> bcs;

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
        if (type == "zeroGradient")
        {
            neoPatchDict.insert("type", std::string("fixedGradient"));
            neoPatchDict.insert("fixedGradient", type_primitive_t {});
        }
        if (type == "extrapolatedCalculated")
        {
            neoPatchDict.insert("type", std::string("calculated"));
        }
        bcs.push_back(fvcc::VolumeBoundary<type_primitive_t>(uMesh, neoPatchDict, patchi));
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

    type_container_t nfVolField(exec, uMesh, readVolBoundaryConditions(uMesh, volField));

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

    std::vector<fvcc::SurfaceBoundary<type_primitive_t>> bcs;

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
        NeoFOAM::Dictionary neoPatchDict;
        neoPatchDict.insert("type", std::string(type));
        bcs.push_back(fvcc::SurfaceBoundary<type_primitive_t>(uMesh, neoPatchDict, patchi));
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
