// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include <t8.h>
#include "NeoFOAM/mesh/treeAMR/treeAMRMesh.hpp"
#include "NeoFOAM/mesh/treeAMR/treeAMRMeshModifier.hpp"

#include "FoamAdapter/treeAMR/t8_MeshModifier.hpp"
#include <memory>

template<typename T>
void printField(const std::vector<T>& field, const std::string& name)
{
    std::cout << name << std::endl;
    for (auto f : field)
    {
        std::cout << f << std::endl;
    }
};

int main(int argc, char** argv)
{
    int mpiret;
    /* The prefix for our output files. */
    // const char prefix[BUFSIZ] = "t8_step1_tetcube";
    // t8_locidx_t local_num_trees;
    // t8_gloidx_t global_num_trees;

    /* Initialize MPI. This has to happen before we initialize sc or t8code. */
    mpiret = sc_MPI_Init(&argc, &argv);
    /* Error check the MPI return value. */
    SC_CHECK_MPI(mpiret);

    sc_MPI_Comm comm = sc_MPI_COMM_WORLD;


    /* Initialize the sc library, has to happen before we initialize t8code. */
    sc_init(sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
    /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels.
     */
    t8_init(SC_LP_PRODUCTION);

    Kokkos::initialize(argc, argv);

    {
        // std::unique_ptr<NeoFOAM::treeAMRMeshModifier> t8mod =
        // std::make_unique<t8_MeshModifier>(comm,5,8);
        NeoFOAM::treeAMRMesh tMesh(std::make_unique<t8_MeshModifier>(comm, 1, 2));

        tMesh.refine();

        tMesh.write();

        std::cout << "nElements: " << tMesh.nElements() << std::endl;

        const NeoFOAM::treeAMRMeshModifier& mod = tMesh.meshModifier();
        // std::cout << "V: " << tMesh.V() << std::endl;
        // print the volume
        printField(mod.V(), "V");
        // print the level
        printField(mod.level(), "level");

        for (int i = 0; i < mod.Cx().size(); i++)
        {
            std::cout << "i: " << i << " Cx: " << mod.Cx()[i] << " Cy: " << mod.Cy()[i]
                      << std::endl;
        }

        // print Ax
        printField(mod.Ax(), "Ax");
        // print Ay
        printField(mod.Ay(), "Ay");
        // print Az
        printField(mod.Az(), "Az");

        // print Cfx
        for (int i = 0; i < mod.Cfx().size(); i++)
        {
            std::cout << "i: " << i << " Cfx: " << mod.Cfx()[i] << " Cfy: " << mod.Cfy()[i]
                      << std::endl;
        }

        printField(mod.owner(), "owner");
        printField(mod.neighbour(), "neighbour");

        printField(mod.level(), "level");
    }

    Kokkos::finalize();

    sc_finalize();

    mpiret = sc_MPI_Finalize();
    SC_CHECK_MPI(mpiret);


    return 0;
}
