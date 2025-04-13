// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2025 NeoFOAM authors

Info << "Reading transportProperties\n" << endl;

Foam::IOdictionary transportProperties(Foam::IOobject(
    "transportProperties",
    runTime.constant(),
    mesh,
    Foam::IOobject::MUST_READ_IF_MODIFIED,
    Foam::IOobject::NO_WRITE
));

Foam::dimensionedScalar nu("nu", Foam::dimViscosity, transportProperties);

Info << "Reading field p\n" << endl;
Foam::volScalarField
    p(Foam::IOobject(
          "p",
          runTime.timeName(),
          mesh,
          Foam::IOobject::MUST_READ,
          Foam::IOobject::AUTO_WRITE
      ),
      mesh);


Info << "Reading field U\n" << endl;
Foam::volVectorField
    U(Foam::IOobject(
          "U",
          runTime.timeName(),
          mesh,
          Foam::IOobject::MUST_READ,
          Foam::IOobject::AUTO_WRITE
      ),
      mesh);


Info << "Reading/calculating face flux field phi\n" << endl;

Foam::surfaceScalarField
    phi(Foam::IOobject(
            "phi",
            runTime.timeName(),
            mesh,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::AUTO_WRITE
        ),
        fvc::flux(U));


Foam::label pRefCell = 0;
Foam::scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
