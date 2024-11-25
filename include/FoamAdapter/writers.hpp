// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "fvMesh.H"
#include "volFields.H"

#include "NeoFOAM/fields/field.hpp"
#include "NeoFOAM/core/error.hpp"

#include "FoamAdapter/conversion/convert.hpp"

namespace Foam
{

namespace detail
{
template<class DestField, class SrcField>
void copy_impl(DestField& dest, const SrcField src)
{
    NF_ASSERT_EQUAL(dest.size(), src.size());
    auto src_host = src.copyToHost();
    auto src_span = src_host.span();
    for (int i = 0; i < dest.size(); i++)
    {
        dest[i] = convert(src_span[i]);
    }
}
}

void write(NeoFOAM::scalarField& sf, const Foam::fvMesh& mesh, const std::string fieldName)
{
    Foam::volScalarField* field = mesh.getObjectPtr<Foam::volScalarField>(fieldName);
    if (field)
    {
        detail::copy_impl(field->ref(), sf);
        field->write();
    }
    else
    {
        Foam::volScalarField foamField(
            Foam::IOobject(
                fieldName,
                mesh.time().timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar(Foam::dimless, 0)
        );
        detail::copy_impl(foamField.ref(), sf);
        foamField.write();
    }
}

void write(NeoFOAM::vectorField& sf, const Foam::fvMesh& mesh, const std::string fieldName)
{
    Foam::volVectorField* field = mesh.getObjectPtr<Foam::volVectorField>(fieldName);
    if (field)
    {
        // field is already present and needs to be updated
        detail::copy_impl(field->ref(), sf);
        field->write();
    }
    else
    {
        Foam::volVectorField foamField(
            Foam::IOobject(
                fieldName,
                mesh.time().timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedVector(Foam::dimless, Foam::Zero)
        );
        detail::copy_impl(foamField.ref(), sf);
        foamField.write();
    }
}

} // namespace Foam
