#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include "VoxelSpaceBase.hpp"

namespace ecell4
{

class LatticeSpace : public VoxelSpaceBase
{
public:

    typedef VoxelSpaceBase base_type;

public:

    LatticeSpace(const Real3& edge_lengths,
                 const Real& voxel_radius,
                 const bool is_periodic);
    virtual ~LatticeSpace();

    void set_lattice_properties(const bool is_periodic);
    void reset(const Real3& edge_lengths,
               const Real& voxel_radius,
               const bool is_periodic);

    const Integer col_size() const;
    const Integer row_size() const;
    const Integer layer_size() const;

    Integer3 coordinate2global(const coordinate_type& coord) const;
    coordinate_type global2coordinate(const Integer3& global) const;

    Real3 global2position(const Integer3& global) const;
    Integer3 position2global(const Real3& pos) const;

    coordinate_type periodic_transpose(const coordinate_type& coord) const;

    bool is_in_range(const coordinate_type& coord) const;
    bool is_inside(const coordinate_type& coord) const;

    /*
     * for LatticeSpaceBase
     */

    Real3 coordinate2position(const coordinate_type& coord) const;
    coordinate_type position2coordinate(const Real3& pos) const;

    Integer num_neighbors(const coordinate_type& coord) const;
    coordinate_type get_neighbor(const coordinate_type& coord, const Integer& nrand) const;
    // virtual coordinate_type
    //     get_neighbor_boundary(const coordinate_type& coord, const Integer& nrand) const = 0;

    Integer size() const;
    Integer3 shape() const;

    // virtual bool on_structure(const Voxel& v) = 0;

    coordinate_type inner2coordinate(const coordinate_type inner) const;
    Integer inner_size() const;

    /*
     * VoxelSpaceBaseTraits
     */

    /*
     * ParticleSpaceTraits
     */
    const Real3& edge_lengths() const;
    Real3 actual_lengths() const;

protected:

    Real3 edge_lengths_;
    Real HCP_L, HCP_X, HCP_Y;
    Integer row_size_, layer_size_, col_size_;

};

} // ecell4

#endif /* __ECELL4_LATTICE_SPACE_HPP */
