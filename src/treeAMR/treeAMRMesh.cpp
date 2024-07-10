// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "treeAMRMesh.hpp"

#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_forest/t8_forest_general.h>            /* forest definition and basic interface. */
#include <t8_forest/t8_forest_io.h>                 /* save forest */
#include <t8_forest/t8_forest_geometrical.h>        /* geometrical information of the forest */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <t8_vec.h>

int t8_step3_adapt_callback(
    t8_forest_t forest,
    t8_forest_t forest_from,
    t8_locidx_t which_tree,
    t8_locidx_t lelement_id,
    t8_eclass_scheme_c* ts,
    const int is_family,
    const int num_elements,
    t8_element_t* elements[]
)
{
    /* Our adaptation criterion is to look at the midpoint coordinates of the current element and if
     * they are inside a sphere around a given midpoint we refine, if they are outside, we coarsen.
     */
    double centroid[3]; /* Will hold the element midpoint. */
    /* In t8_step3_adapt_forest we pass a t8_step3_adapt_data pointer as user data to the
     * t8_forest_new_adapt function. This pointer is stored as the used data of the new forest
     * and we can now access it with t8_forest_get_user_data (forest). */
    const struct t8_step3_adapt_data* adapt_data =
        (const struct t8_step3_adapt_data*)t8_forest_get_user_data(forest);
    double dist; /* Will store the distance of the element's midpoint and the sphere midpoint. */

    /* You can use T8_ASSERT for assertions that are active in debug mode (when configured with
     * --enable-debug). If the condition is not true, then the code will abort. In this case, we
     * want to make sure that we actually did set a user pointer to forest and thus did not get the
     * NULL pointer from t8_forest_get_user_data.
     */
    T8_ASSERT(adapt_data != NULL);

    /* Compute the element's centroid coordinates. */
    t8_forest_element_centroid(forest_from, which_tree, elements[0], centroid);

    /* Compute the distance to our sphere midpoint. */
    dist = t8_vec_dist(centroid, adapt_data->midpoint);
    if (dist < adapt_data->refine_if_inside_radius)
    {
        /* Refine this element. */
        return 1;
    }
    else if (is_family && dist > adapt_data->coarsen_if_outside_radius)
    {
        /* Coarsen this family. Note that we check for is_family before, since returning < 0
         * if we do not have a family as input is illegal. */
        return -1;
    }
    /* Do not change this element. */
    return 0;
}

// Implement the member functions of the treeAMRMesh class here

treeAMRMesh::treeAMRMesh(sc_MPI_Comm comm, int initialLevel_, int maxLevel_)
    : initialLevel_(initialLevel_), maxLevel_(maxLevel_), nElements_(0)
{
    /* Initialize an adapted forest with periodic boundaries. */
    cmesh_ = t8_cmesh_new_hypercube(T8_ECLASS_QUAD, comm, 0, 0, 0);

    scheme = t8_scheme_new_default_cxx(); // deleting the pointer will causes sc abort to fail
    forest_ = t8_forest_new_uniform(cmesh_, scheme, initialLevel_, 0, comm);
    nElements_ = t8_forest_get_local_num_elements(forest_);
}

treeAMRMesh::~treeAMRMesh()
{
    // Destructor code goes here
    t8_cmesh_destroy(&cmesh_);
    t8_forest_unref(&forest_);
}

bool treeAMRMesh::refine()
{
    t8_forest_t forest_adapt;
    struct t8_step3_adapt_data adapt_data = {
        {0.0, 0.0, 0}, /* Midpoints of the sphere. */
        0.4,           /* Refine if inside this radius. */
        2.0            /* Coarsen if outside this radius. */
    };

    /* Check that forest is a committed, that is valid and usable, forest. */
    T8_ASSERT(t8_forest_is_committed(forest));

    t8_forest_t tmp_forest;

    t8_forest_init(&tmp_forest);
    t8_forest_set_user_data(tmp_forest, &adapt_data);
    t8_forest_set_adapt(tmp_forest, forest_, t8_step3_adapt_callback, 0);
    t8_forest_set_partition(tmp_forest, NULL, 0);
    t8_forest_set_balance(tmp_forest, NULL, 0);
    t8_forest_set_ghost(tmp_forest, 1, T8_GHOST_FACES);
    t8_forest_commit(tmp_forest);

    forest_ = tmp_forest;
    nElements_ = t8_forest_get_local_num_elements(forest_);

    return true;
}

bool treeAMRMesh::coarsen()
{
    // Coarsening code goes here
    return true;
}

bool treeAMRMesh::balance()
{
    // Balancing code goes here
    return true;
}

void treeAMRMesh::write()
{
    // Write code goes here
    t8_cmesh_vtk_write_file(cmesh_, "prefix", 1.0);
    t8_forest_write_vtk(forest_, "forest");
}

int32_t treeAMRMesh::nElements() const { return nElements_; }
