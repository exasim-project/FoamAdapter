// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors


#include "fvCFD.H"
// #include "profiling.H"
// #include "benchmark/benchmark.h"
// #include <nanobench.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// void BM_TEqnTemporalOperator(volScalarField& T) {
//     fvScalarMatrix TEqnTemporalOperator(fvm::ddt(T));
// }

// void BM_TEqnDiffusionOperator(volScalarField& Gamma,volScalarField& T) {
//     fvScalarMatrix TEqnDiffusionOperator(fvm::laplacian(Gamma, T));
// }

// void BM_TEqnConvectionOperator(surfaceScalarField& phi, volScalarField& T) {
//     fvScalarMatrix TEqnConvectionOperator(fvm::div(phi, T));
// }

// void BM_EnergyEquation(volScalarField& rho,surfaceScalarField& phi, volScalarField&
// Gamma,volScalarField& T) {
//     fvScalarMatrix EnergyEquation(fvm::ddt(rho, T) + fvm::div(phi, T) - fvm::laplacian(Gamma,
//     T));
// }

// void gen(std::string const& typeName, char const* mustacheTemplate,

//     ankerl::nanobench::Bench const& bench) {


//     std::ofstream templateOut("mustache.template." + typeName);

//     templateOut << mustacheTemplate;


//     std::ofstream renderOut("mustache.render." + typeName);

//     ankerl::nanobench::render(mustacheTemplate, bench, renderOut);

// }

int main(int argc, char* argv[])
{
#include "addProfilingOption.H"
#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
    // ankerl::nanobench::Bench bench;
    // bench.title("matrixAssembly");

    // bench.run("BM_TEqnTemporalOperator", [&] {
    //     BM_TEqnTemporalOperator(T);
    // });

    // bench.run("BM_TEqnDiffusionOperator", [&] {
    //     BM_TEqnDiffusionOperator(Gamma,T);
    // });

    // bench.run("BM_TEqnConvectionOperator", [&] {
    //     BM_TEqnConvectionOperator(phi,T);
    // });

    //     bench.run("BM_EnergyEquation", [&] {
    //     BM_EnergyEquation(rho,phi,Gamma,T);
    // });

    // gen("json", ankerl::nanobench::templates::json(), bench);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
