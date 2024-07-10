// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/mesh/treeAMR/treeAMRMeshModifier.hpp"
#include <cstdint>

#include <t8.h>                              /* General t8code header, always include this. */
#include <t8_cmesh.h>                        /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>      /* A collection of exemplary cmeshes */
#include <t8_cmesh_vtk_writer.h>             /* cmesh-writer interface. */
#include <t8_forest/t8_forest_general.h>     /* forest definition and basic interface. */
#include <t8_forest/t8_forest_geometrical.h> /* geometrical information of the forest */
#include <t8_forest/t8_forest_io.h>          /* save forest */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <t8_vec.h>                                 /* Basic operations on 3D vectors. */

struct t8_step3_adapt_data
{
    double midpoint[3];               /* The midpoint of our sphere. */
    double refine_if_inside_radius;   /* if an element's center is smaller than this
                                         value, we refine the element. */
    double coarsen_if_outside_radius; /* if an element's center is larger this
                                         value, we coarsen its family. */
};

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
    /* Our adaptation criterion is to look at the midpoint coordinates of the
     * current element and if they are inside a sphere around a given midpoint we
     * refine, if they are outside, we coarsen. */
    double centroid[3]; /* Will hold the element midpoint. */
    /* In t8_step3_adapt_forest we pass a t8_step3_adapt_data pointer as user data
     * to the t8_forest_new_adapt function. This pointer is stored as the used
     * data of the new forest and we can now access it with
     * t8_forest_get_user_data (forest). */
    const struct t8_step3_adapt_data* adapt_data =
        (const struct t8_step3_adapt_data*)t8_forest_get_user_data(forest);
    double dist; /* Will store the distance of the element's midpoint and the
                    sphere midpoint. */

    /* You can use T8_ASSERT for assertions that are active in debug mode (when
     * configured with --enable-debug). If the condition is not true, then the
     * code will abort. In this case, we want to make sure that we actually did
     * set a user pointer to forest and thus did not get the NULL pointer from
     * t8_forest_get_user_data.
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
        /* Coarsen this family. Note that we check for is_family before, since
         * returning < 0 if we do not have a family as input is illegal. */
        return -1;
    }
    /* Do not change this element. */
    return 0;
};

/**
 * @class t8_MeshModifier
 * @brief Represents a mesh modifier for tree-based adaptive mesh refinement
 * (AMR).
 *
 * The t8_MeshModifier class inherits from the treeAMRMeshModifier class and
 * provides additional functionality for refining, coarsening, balancing, and
 * writing the mesh. It also provides information about the number of elements
 * in the mesh.
 */
class t8_MeshModifier : public NeoFOAM::treeAMRMeshModifier
{
public:

    /**
     * @brief Constructor for treeAMRMesh.
     * @param comm The MPI communicator.
     * @param initialLevel_ The initial level of refinement for the mesh.
     * @param maxLevel_ The maximum level of refinement for the mesh.
     */
    t8_MeshModifier(sc_MPI_Comm comm, int32_t initialLevel_, int32_t maxLevel_)

        : treeAMRMeshModifier(initialLevel_, maxLevel_)
    {
        /* Initialize an adapted forest with periodic boundaries. */
        cmesh_ = t8_cmesh_new_hypercube(T8_ECLASS_QUAD, comm, 0, 0, 0);

        scheme = t8_scheme_new_default_cxx(); // deleting the pointer will causes sc
                                              // abort to fail
        forest_ = t8_forest_new_uniform(cmesh_, scheme, initialLevel_, 0, comm);
        nElements_ = t8_forest_get_local_num_elements(forest_);
    };

    /**
     * @brief Destructor for t8_MeshModifier.
     */
    virtual ~t8_MeshModifier()
    {
        // Destructor code goes here
        t8_cmesh_destroy(&cmesh_);
        t8_forest_unref(&forest_);
    };

    /**
     * @brief Refines the mesh.
     * @return True if the mesh was successfully refined, false otherwise.
     */
    virtual bool refine()
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

        computeFaceAddressing();

        return true;
    };

    /**
     * @brief Coarsens the mesh.
     * @return True if the mesh was successfully coarsened, false otherwise.
     */
    virtual bool coarsen()
    {
        // Coarsening code goes here
        return true;
    };

    /**
     * @brief Balances the mesh.
     * @return True if the mesh was successfully balanced, false otherwise.
     */
    virtual bool balance()
    {
        // Balancing code goes here
        return true;
    };

    /**
     * @brief Writes the mesh to a file.
     */
    virtual void write()
    {
        // Write code goes here
        t8_cmesh_vtk_write_file(cmesh_, "cmesh", 1.0);
        t8_forest_write_vtk(forest_, "forest");
    };

    /**
     * @brief Gets the number of elements in the mesh.
     * @return The number of elements in the mesh.
     */
    int32_t nElements() const { return nElements_; };

private:

    int32_t initialLevel_; ///< The initial level of refinement for the mesh.
    int32_t maxLevel_;     ///< The maximum level of refinement for the mesh.
    int32_t nElements_;    ///< The number of elements in the mesh.

    t8_cmesh_t cmesh_;       ///< The underlying mesh data structure.
    t8_forest_t forest_;     ///< The underlying forest data structure.
    t8_scheme_cxx_t* scheme; ///< The scheme used for mesh operations.

    void computeFaceAddressing()
    {
        std::array<double, 3> elem_midpoint = {0.0, 0.0, 0.0};
        std::array<double, 3> face_normal = {0.0, 0.0, 0.0};
        std::array<double, 3> face_centroid = {0.0, 0.0, 0.0};
        std::array<double, 3> dist_face_to_center = {0.0, 0.0, 0.0};
        double V = 0.0;

        // loop over the elements of forests
        int32_t lclIdx = 0;
        for (int32_t itree = 0; itree < t8_forest_get_num_local_trees(forest_); itree++)
        {
            t8_eclass_t tree_class = t8_forest_get_tree_class(forest_, itree);
            t8_eclass_scheme_c* eclass_scheme = t8_forest_get_eclass_scheme(forest_, tree_class);
            for (int32_t ielement = 0; ielement < t8_forest_get_tree_num_elements(forest_, itree);
                 ielement++)
            {

                const t8_element_t* element =
                    t8_forest_get_element_in_tree(forest_, itree, ielement);
                t8_forest_element_centroid(forest_, itree, element, elem_midpoint.data());
                V = t8_forest_element_volume(forest_, itree, element);
                int level = eclass_scheme->t8_element_level(element);
                Cx_.push_back(elem_midpoint[0]);
                Cy_.push_back(elem_midpoint[1]);
                Cz_.push_back(elem_midpoint[2]);
                V_.push_back(V);
                level_.push_back(level);

                // eclass_scheme->t8_element_boundary_face

                int num_faces = eclass_scheme->t8_element_num_faces(element);
                for (int iface = 0; iface < num_faces; iface++)
                {
                    int num_neighbors;        /**< Number of neighbors for each face */
                    int* dual_faces;          /**< The face indices of the neighbor elements */
                    t8_locidx_t* neighids;    /**< Indices of the neighbor elements */
                    t8_element_t** neighbors; /*< Neighboring elements. */
                    t8_eclass_scheme_c* neigh_scheme; /*< Neighboring elements scheme. */

                    /* Collect all neighbors at the current face. */
                    t8_forest_leaf_face_neighbors(
                        forest_,
                        itree,
                        element,
                        &neighbors,
                        iface,
                        &dual_faces,
                        &num_neighbors,
                        &neighids,
                        &neigh_scheme,
                        1
                    );
                    std::cout << "iface: " << iface << std::endl;
                    if (num_neighbors > 0)
                    {

                        std::cout << "num_neighbors: " << num_neighbors << std::endl;

                        for (int ineigh = 0; ineigh < num_neighbors; ineigh++)
                        {
                            std::cout << "dual_faces: " << dual_faces[ineigh] << std::endl;
                        }

                        t8_forest_element_face_normal(
                            forest_, itree, element, iface, face_normal.data()
                        );
                        t8_forest_element_face_centroid(
                            forest_, itree, element, iface, face_centroid.data()
                        );
                        // t8_element_level()
                        // does the normal point inwards or outwards?
                        dist_face_to_center[0] = elem_midpoint[0] - face_centroid[0];
                        dist_face_to_center[1] = elem_midpoint[1] - face_centroid[1];
                        dist_face_to_center[2] = elem_midpoint[2] - face_centroid[2];
                        // can have multiple neibouts
                        // TODO account for multiple neighbors
                        Cfx_.push_back(face_centroid[0]);
                        Cfy_.push_back(face_centroid[1]);
                        Cfz_.push_back(face_centroid[2]);
                        std::cout << "face_normal: " << face_normal[0] << " " << face_normal[1]
                                  << " " << face_normal[2] << std::endl;

                        double orient_x = dist_face_to_center[0] * face_normal[0];
                        double orient_y = dist_face_to_center[1] * face_normal[1];
                        double orient_z = dist_face_to_center[2] * face_normal[2];
                        // get index of absolute max value

                        // Find the iterator of the largest absolute value
                        auto max_it = std::max_element(
                            face_normal.begin(),
                            face_normal.end(),
                            [](double a, double b) { return std::abs(a) < std::abs(b); }
                        );

                        // Calculate the index
                        int orient = std::distance(face_normal.begin(), max_it);

                        // only face with positive normal are stored
                        // otherwise we would store the face twice
                        if (face_normal[orient] > 0.0)
                        {
                            owner_.push_back(lclIdx);
                            neighbour_.push_back(neighids[0]);
                            switch (orient)
                            {
                            case 0:
                                Ax_.push_back(
                                    t8_forest_element_face_area(forest_, itree, element, iface)
                                );
                                break;
                            case 1:
                                Ay_.push_back(
                                    t8_forest_element_face_area(forest_, itree, element, iface)
                                );
                                break;
                            case 2:
                                Az_.push_back(
                                    t8_forest_element_face_area(forest_, itree, element, iface)
                                );
                                break;
                            }
                        }
                        // else if(orient_y < 0.0)
                        // {
                        //     Ay_.push_back(t8_forest_element_face_area(forest_, itree,
                        //     element, iface));
                        // }
                        // else if(orient_z < 0.0)
                        // {
                        //     Az_.push_back(t8_forest_element_face_area(forest_, itree,
                        //     element, iface));
                        // }
                    }

                    /* Retrieve the `height` of the face neighbor. Account for two
                    neighbors in case of a non-conforming interface by computing the
                    average. */
                    // double height = 0.0;
                    // if (num_neighbors > 0) {
                    // for (int ineigh = 0; ineigh < num_neighbors; ineigh++) {
                    //     height = height + element_data[neighids[ineigh]].height;
                    // }
                    // height = height / num_neighbors;
                    // }

                    // /* Fill in the neighbor information of the 3x3 stencil. */
                    // switch (iface) {
                    // case 0:  // NORTH
                    // stencil[0][1] = height;
                    // dx[0] = element_data[neighids[0]].dx;
                    // break;
                    // case 1:  // SOUTH
                    // stencil[2][1] = height;
                    // dx[2] = element_data[neighids[0]].dx;
                    // break;
                    // case 2:  // WEST
                    // stencil[1][0] = height;
                    // dy[0] = element_data[neighids[0]].dy;
                    // break;
                    // case 3:  // EAST
                    // stencil[1][2] = height;
                    // dy[2] = element_data[neighids[0]].dy;
                    // break;
                    // }

                    /* Free allocated memory. */
                    T8_FREE(neighbors);
                    T8_FREE(dual_faces);
                    T8_FREE(neighids);
                }
                // increment the local index over multiple trees
                lclIdx++;
            }
        }
    };
};
