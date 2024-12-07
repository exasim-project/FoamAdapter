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
auto readVolBoundaryConditions(const NeoFOAM::UnstructuredMesh& nfMesh, const FoamType& ofVolField)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;

    // get boundary as dictionary
    OStringStream os;
    ofVolField.boundaryField().writeEntries(os);
    IStringStream is(os.str());
    dictionary bDict(is);

    std::map<std::string, std::function<void(NeoFOAM::Dictionary&)>> patchInserter {
        {"fixedGradient", [](auto& dict) { dict.insert("type", std::string("fixedGradient")); }},
        {"zeroGradient",
         [&](auto& dict)
         {
             dict.insert("type", std::string("fixedGradient"));
             dict.insert("fixedGradient", type_primitive_t {});
         }},
        {"fixedValue",
         [](auto& dict)
         {
             dict.insert("type", std::string("fixedValue"));
             dict.insert("fixedValue", type_primitive_t {});
         }},
        {"calculated", [](auto& dict) { dict.insert("type", std::string("calculated")); }},
        {"extrapolatedCalculated",
         [](auto& dict) { dict.insert("type", std::string("calculated")); }},
        {"empty", [](auto& dict) { dict.insert("type", std::string("empty")); }}
    };

    int patchi = 0;
    std::vector<fvcc::VolumeBoundary<type_primitive_t>> bcs;
    for (const auto& bName : bDict.toc())
    {
        dictionary patchDict = bDict.subDict(bName);
        NeoFOAM::Dictionary neoPatchDict;
        patchInserter[patchDict.get<word>("type")](neoPatchDict);
        bcs.emplace_back(nfMesh, neoPatchDict, patchi);
        patchi++;
    }
    return bcs;
}

template<typename FoamType>
auto constructFrom(
    const NeoFOAM::Executor exec, const NeoFOAM::UnstructuredMesh& nfMesh, const FoamType& in
)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;

    type_container_t out(exec, in.name(), nfMesh, readVolBoundaryConditions(nfMesh, in));

    out.internalField() = fromFoamField(exec, in.primitiveField());
    out.correctBoundaryConditions();

    return out;
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
    OStringStream os;
    surfaceField.boundaryField().writeEntries(os);
    IStringStream is(os.str());
    dictionary bDict(is);
    int patchi = 0;

    std::map<std::string, std::function<void(NeoFOAM::Dictionary&)>> patchInserter {
        {"fixedGradient", [](auto& dict) { dict.insert("type", std::string("fixedGradient")); }},
        {"zeroGradient",
         [&](auto& dict)
         {
             dict.insert("type", std::string("fixedGradient"));
             dict.insert("fixedGradient", type_primitive_t {});
         }},
        {"fixedValue",
         [](auto& dict)
         {
             dict.insert("type", std::string("fixedValue"));
             dict.insert("fixedValue", type_primitive_t {});
         }},
        {"calculated", [](auto& dict) { dict.insert("type", std::string("calculated")); }},
        {"empty", [](auto& dict) { dict.insert("type", std::string("empty")); }}
    };

    for (const auto& bName : bDict.toc())
    {
        dictionary patchDict = bDict.subDict(bName);
        NeoFOAM::Dictionary neoPatchDict;
        patchInserter[patchDict.get<word>("type")](neoPatchDict);
        bcs.push_back(fvcc::SurfaceBoundary<type_primitive_t>(uMesh, neoPatchDict, patchi));
        patchi++;
    }
    return bcs;
}

template<typename FoamType>
auto constructSurfaceField(
    const NeoFOAM::Executor exec, const NeoFOAM::UnstructuredMesh& nfMesh, const FoamType& in
)
{
    using type_container_t = typename type_map<FoamType>::container_type;
    using type_primitive_t = typename type_map<FoamType>::mapped_type;
    using foam_primitive_t = typename FoamType::cmptType;

    type_container_t out(exec, in.name(),nfMesh, std::move(readSurfaceBoundaryConditions(nfMesh, in)));

    Field<foam_primitive_t> flattenedField(out.internalField().size());
    size_t nInternal = nfMesh.nInternalFaces();

    forAll(in, facei)
    {
        flattenedField[facei] = convert(in[facei]);
    }

    label idx = nInternal;
    forAll(in.boundaryField(), patchi)
    {
        const fvsPatchField<foam_primitive_t>& pin = in.boundaryField()[patchi];

        forAll(pin, facei)
        {
            flattenedField[idx] = pin[facei];
            idx++;
        }
    }
    assert(idx == flattenedField.size());

    out.internalField() = fromFoamField(exec, flattenedField);
    out.correctBoundaryConditions();

    return out;
}

}; // namespace Foam
