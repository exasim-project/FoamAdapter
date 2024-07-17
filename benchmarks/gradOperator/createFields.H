// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

Foam::Info << "Reading field T\n" << Foam::endl;

Foam::volScalarField
    T(Foam::IOobject(
          "T", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
      ),
      mesh);


Foam::Info << "Reading diffusivity Gamma\n" << Foam::endl;

Foam::volScalarField Gamma(
    Foam::IOobject(
        "Gamma", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
    ),
    mesh
);


Foam::Info << "Reading diffusivity Gamma\n" << Foam::endl;

Foam::volScalarField
    rho(Foam::IOobject(
            "rho", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
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
    linearInterpolate(rho* U) & mesh.Sf()
);

Foam::volScalarField
    mu(Foam::IOobject(
           "mu", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::AUTO_WRITE
       ),
       mesh);
