
#ifndef TREEAMRMESH_HPP
#define TREEAMRMESH_HPP

#include <cstdint>

#include <t8.h>
#include <t8_cmesh.h>
#include <t8_forest/t8_forest_general.h> /* forest definition and general interface. */
#include <t8_forest/t8_forest_io.h>                 /* forest io interface. */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */

struct t8_step3_adapt_data {
  double midpoint[3];             /* The midpoint of our sphere. */
  double refine_if_inside_radius; /* if an element's center is smaller than this
                                     value, we refine the element. */
  double coarsen_if_outside_radius; /* if an element's center is larger this
                                       value, we coarsen its family. */
};

int t8_step3_adapt_callback(t8_forest_t forest, t8_forest_t forest_from,
                            t8_locidx_t which_tree, t8_locidx_t lelement_id,
                            t8_eclass_scheme_c *ts, const int is_family,
                            const int num_elements, t8_element_t *elements[]);

/**
 * @class treeAMRMesh
 * @brief Represents a mesh for tree-based adaptive mesh refinement (AMR).
 *
 * The treeAMRMesh class provides functionality for refining, coarsening,
 * balancing, and writing the mesh. It also provides information about the
 * number of elements in the mesh.
 */
class treeAMRMesh {
public:
  /**
   * @brief Constructor for treeAMRMesh.
   * @param comm The MPI communicator.
   * @param initialLevel_ The initial level of refinement for the mesh.
   * @param maxLevel_ The maximum level of refinement for the mesh.
   */
  treeAMRMesh(sc_MPI_Comm comm, int32_t initialLevel_, int32_t maxLevel_);

  /**
   * @brief Destructor for treeAMRMesh.
   */
  ~treeAMRMesh();

  /**
   * @brief Refines the mesh.
   * @return True if the mesh was successfully refined, false otherwise.
   */
  bool refine();

  /**
   * @brief Coarsens the mesh.
   * @return True if the mesh was successfully coarsened, false otherwise.
   */
  bool coarsen();

  /**
   * @brief Balances the mesh.
   * @return True if the mesh was successfully balanced, false otherwise.
   */
  bool balance();

  /**
   * @brief Writes the mesh to a file.
   */
  void write();

  /**
   * @brief Gets the number of elements in the mesh.
   * @return The number of elements in the mesh.
   */
  int32_t nElements() const;

private:
  int32_t initialLevel_; ///< The initial level of refinement for the mesh.
  int32_t maxLevel_;     ///< The maximum level of refinement for the mesh.
  int32_t nElements_;    ///< The number of elements in the mesh.

  t8_cmesh_t cmesh_;       ///< The underlying mesh data structure.
  t8_forest_t forest_;     ///< The underlying forest data structure.
  t8_scheme_cxx_t *scheme; ///< The scheme used for mesh operations.
};

#endif // TREEAMRMESH_HPP
