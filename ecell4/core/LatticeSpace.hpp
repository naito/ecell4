#ifndef ECELL4_LATTICE_SPACE_HPP
#define ECELL4_LATTICE_SPACE_HPP

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

    Real unit_voxel_volume() const;

    /*
     * ParticleSpaceTraits
     */
    const Real3& edge_lengths() const;
    Real3 actual_lengths() const;

    /*
     * static members
     */
    static inline Real calculate_voxel_volume(const Real voxel_radius)
    {
        return 4.0 * sqrt(2.0) * pow(voxel_radius, 3);
    }

    static inline Real3 calculate_hcp_lengths(const Real voxel_radius)
    {
        return Real3(
                voxel_radius / sqrt(3.0),
                voxel_radius * sqrt(8.0 / 3.0),
                voxel_radius * sqrt(3.0));
    }

    static inline
    Integer3
    calculate_shape(const Real3& edge_lengths, const Real& voxel_radius, const bool is_periodic)
    {
        const Real3 hcpLXY = calculate_hcp_lengths(voxel_radius);
        const Real lengthX = edge_lengths[0];
        const Real lengthY = edge_lengths[1];
        const Real lengthZ = edge_lengths[2];

        Integer col_size = (Integer)rint(lengthX / hcpLXY[1]) + 1;
        Integer layer_size = (Integer)rint(lengthY / hcpLXY[2]) + 1;
        Integer row_size = (Integer)rint((lengthZ / 2) / voxel_radius) + 1;

        if (is_periodic)
        {
            // The number of voxels in each axis must be even for a periodic boundary.
            col_size = (col_size % 2 == 0 ? col_size : col_size + 1);
            layer_size = (layer_size % 2 == 0 ? layer_size : layer_size + 1);
            row_size = (row_size % 2 == 0 ? row_size : row_size + 1);
        }

        return Integer3(col_size, row_size, layer_size);
    }

    static inline
    Real
    calculate_volume(const Real3& edge_lengths, const Real& voxel_radius, const bool is_periodic)
    {
        const Integer3 shape = calculate_shape(edge_lengths, voxel_radius, is_periodic);
        return static_cast<Real>(shape[0] * shape[1] * shape[2])
                * calculate_voxel_volume(voxel_radius);
    }

protected:

    Real3 edge_lengths_;
    Real HCP_L, HCP_X, HCP_Y;
    Integer row_size_, layer_size_, col_size_;

};

} // ecell4

#endif /* ECELL4_LATTICE_SPACE_HPP */
