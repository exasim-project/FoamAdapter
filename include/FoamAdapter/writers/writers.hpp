// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "FoamAdapter/conversion/convert.hpp"
#include "NeoFOAM/fields/FieldTypeDefs.hpp"
#include "fvMesh.H"
#include "volFields.H"

namespace Foam
{

void write(NeoFOAM::scalarField& sf, const Foam::fvMesh& mesh, const std::string fieldName)
{
    Foam::volScalarField* field = mesh.getObjectPtr<Foam::volScalarField>(fieldName);
    if (field)
    {
        // field is already present and needs to be updated
        auto sf_host = sf.copyToHost().field();
        Foam::scalarField& field_ref = field->ref();
        for (int i = 0; i < field_ref.size(); i++)
        {
            field_ref[i] = sf_host[i];
        }
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
        auto sf_host = sf.copyToHost().field();
        Foam::scalarField& field_ref = foamField.ref();
        for (int i = 0; i < field_ref.size(); i++)
        {
            field_ref[i] = sf_host[i];
        }
        foamField.write();
    }
}

void write(NeoFOAM::vectorField& sf, const Foam::fvMesh& mesh, const std::string fieldName)
{
    Foam::volVectorField* field = mesh.getObjectPtr<Foam::volVectorField>(fieldName);
    if (field)
    {
        // field is already present and needs to be updated
        auto sf_host = sf.copyToHost().field();
        Foam::vectorField& field_ref = field->ref();
        for (int i = 0; i < field_ref.size(); i++)
        {
            field_ref[i] = convert(sf_host[i]);
        }
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
        auto sf_host = sf.copyToHost().field();
        Foam::vectorField& field_ref = foamField.ref();
        for (int i = 0; i < field_ref.size(); i++)
        {
            field_ref[i] = convert(sf_host[i]);
        }
        foamField.write();
    }
}

} // namespace Foam
