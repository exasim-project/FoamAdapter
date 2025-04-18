// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "fvMesh.H"
#include "volFields.H"

#include "NeoN/fields/field.hpp"
#include "NeoN/core/error.hpp"

#include "FoamAdapter/conversion/convert.hpp"

namespace Foam
{

namespace detail
{
template<class DestField, class SrcField>
void copy_impl(DestField& dest, const SrcField src)
{
    NF_ASSERT_EQUAL(dest.size(), src.size());
    auto srcHost = src.copyToHost();
    auto srcView = srcHost.view();
    for (int i = 0; i < dest.size(); i++)
    {
        dest[i] = convert(srcView[i]);
    }
}
}

void write(NeoN::scalarField& sf, const Foam::fvMesh& mesh, const std::string fieldName)
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

void write(NeoN::vectorField& sf, const Foam::fvMesh& mesh, const std::string fieldName)
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
