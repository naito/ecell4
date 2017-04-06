#ifndef __ECELL4_OFFLATTICE_SPACE_HPP
#define __ECELL4_OFFLATTICE_SPACE_HPP

#include "VoxelSpaceBase.hpp"

namespace ecell4
{

class OffLatticeSpace : public VoxelSpaceBase
{
protected:

    typedef VoxelSpaceBase base_type;
    typedef std::vector<VoxelPool*> voxel_container;
    typedef std::vector<std::vector<coordinate_type> > adjoining_container;

public:

    typedef std::pair<coordinate_type, coordinate_type> coordinate_pair_type;
    typedef std::vector<coordinate_pair_type> coordinate_pair_list_type;
    typedef std::vector<Real3> position_container;

public:
    OffLatticeSpace(const Real& voxel_radius);
    ~OffLatticeSpace();

    void reset(const position_container& positions,
               const coordinate_pair_list_type& adjoining_pairs);

    void set_unit_voxel_volume(const Real& volume);

    /*
     * VoxelSpaceBaseTraits
     */
    Real unit_voxel_volume() const;

    identified_voxel get_voxel_at(const coordinate_type& coord) const;
    const VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const;
          VoxelPool* get_voxel_pool_at(const coordinate_type& coord);
    const Particle particle_at(const coordinate_type& coord) const;

    bool update_voxel(const ParticleID& pid, const Voxel& v);
    bool remove_voxel(const ParticleID& pid);
    bool remove_voxel(const coordinate_type& coord);

    bool can_move(const coordinate_type& src,
                          const coordinate_type& dest) const;
    bool move(const coordinate_type& src,
              const coordinate_type& dest,
              const std::size_t candidate=0);
    std::pair<coordinate_type, bool>
        move_to_neighbor(VoxelPool* const& from,
                         VoxelPool* const& loc,
                         coordinate_id_pair_type& info,
                         const Integer nrand);

    /*
     * for LatticeSpaceBase
     */
    Real3 coordinate2position(const coordinate_type& coord) const;
    coordinate_type position2coordinate(const Real3& pos) const;

    Integer num_neighbors(const coordinate_type& coord) const;
    coordinate_type get_neighbor(const coordinate_type& coord, const Integer& nrand) const;
    coordinate_type get_neighbor_boundary(const coordinate_type& coord, const Integer& nrand) const;

    Integer size() const;
    Integer3 shape() const;

    coordinate_type inner2coordinate(const coordinate_type inner) const;
    Integer inner_size() const;

    // Real3 actual_lengths() const;

// #ifdef WITH_HDF5
//     void save_hdf5(H5::Group* root) const;
//     void load_hdf5(const H5::Group& root);
// #endif

protected:

    Integer count_voxels(const boost::shared_ptr<VoxelPool>& vp) const;

    bool is_in_range(const coordinate_type& coord) const;
    coordinate_type get_coord(const ParticleID& pid) const;

protected:

    voxel_container voxels_;
    position_container positions_;
    adjoining_container adjoinings_;
    Real unit_voxel_volume_;
};

boost::shared_ptr<OffLatticeSpace>
create_cubic_offlattice_space(const Real voxel_radius, const Integer3& lattice_size);

boost::shared_ptr<OffLatticeSpace>
create_hcp_offlattice_space(const Real voxel_radius, const Integer3& lattice_size);

} // ecell4

#endif /* __ECELL4_OFFLATTICE_SPACE_HPP */
