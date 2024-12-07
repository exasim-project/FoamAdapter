// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

Foam::Info << "Reading field T\n" << Foam::endl;

Foam::volScalarField
    T(Foam::IOobject(
          "T", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
      ),
      mesh);

Foam::volVectorField
    U(Foam::IOobject(
          "U", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
      ),
      mesh);

Foam::surfaceScalarField phi(
    Foam::IOobject(
        "phi", runTime.timeName(), mesh, Foam::IOobject::READ_IF_PRESENT, Foam::IOobject::AUTO_WRITE
    ),
    Foam::linearInterpolate(U) & mesh.Sf()
);

// Copies of initial U and phi for use when flow is periodic
Foam::surfaceScalarField phi0 = phi;
Foam::volVectorField U0 = U;

