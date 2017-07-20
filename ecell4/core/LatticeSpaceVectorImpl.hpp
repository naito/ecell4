#ifndef __ECELL4_LATTICE_SPACE_VECTOR_IMPL_HPP
#define __ECELL4_LATTICE_SPACE_VECTOR_IMPL_HPP

#include "LatticeSpace.hpp"

namespace ecell4 {

class LatticeSpaceVectorImpl : public LatticeSpace
{
public:

    typedef LatticeSpace base_type;
    typedef std::vector<VoxelPool*> voxel_container;

public:

    LatticeSpaceVectorImpl(
        const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic = true);
    ~LatticeSpaceVectorImpl();

    /*
     * Space APIs
     *
     * using ParticleID, Species and Posision3
     */

    virtual bool remove_voxel(const ParticleID& pid);
    virtual bool remove_voxel(const coordinate_type& coord);

    bool update_structure(const Particle& p);

    /*
     * for Simulator
     *
     * using Species and coordinate_type
     */
    virtual identified_voxel get_voxel_at(const coordinate_type& coord) const;

    // virtual void update_voxel(const Voxel& v);
    virtual bool update_voxel(const ParticleID& pid, const Voxel& v);

    bool add_voxels(const Species species, std::vector<std::pair<ParticleID, coordinate_type> > voxels);

    const Species& find_species(std::string name) const;
    std::vector<coordinate_type> list_coords(const Species& sp) const;
    std::vector<coordinate_type> list_coords_exact(const Species& sp) const;

    virtual const VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const;
    virtual       VoxelPool* get_voxel_pool_at(const coordinate_type& coord);

    // bool update_molecule(coordinate_type coord, const Species& species);
    // bool add_molecule(const Species& sp, coordinate_type coord, const ParticleID& pid);
    virtual bool move(
        const coordinate_type& src, const coordinate_type& dest,
        const std::size_t candidate=0);
    virtual bool can_move(const coordinate_type& src, const coordinate_type& dest) const;

    coordinate_type get_neighbor_boundary(
        const coordinate_type& coord, const Integer& nrand) const
    {
        coordinate_type const dest = get_neighbor(coord, nrand);
        VoxelPool* dest_vp(voxels_.at(dest));
        return (dest_vp != periodic_ ? dest : periodic_transpose(dest));
    }

    inline bool is_periodic() const
    {
        return is_periodic_;
    }

#ifdef WITH_HDF5
    /*
     * HDF5 Save
     */
    void save_hdf5(H5::Group* root) const
    {
        save_lattice_space(*this, root);
    }

    void load_hdf5(const H5::Group& root)
    {
        load_lattice_space(root, this);
    }
#endif

    void reset(const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic)
    {
        base_type::reset(edge_lengths, voxel_radius, is_periodic);

        is_periodic_ = is_periodic;
        initialize_voxels(is_periodic_);
    }

    virtual const Particle particle_at(const coordinate_type& coord) const;

    coordinate_type apply_boundary_(
        const coordinate_type& coord) const
    {
        return periodic_transpose(coord);
    }

protected:

    void initialize_voxels(const bool is_periodic);

    std::pair<coordinate_type, bool> move_(
            coordinate_type from, coordinate_type to,
            const std::size_t candidate=0);
    std::pair<coordinate_type, bool> move_(
            coordinate_id_pair_type& info, coordinate_type to);
    coordinate_type get_coord(const ParticleID& pid) const;

    Integer count_voxels(const boost::shared_ptr<VoxelPool>& vp) const;

protected:

    bool is_periodic_;

    voxel_container voxels_;

    VoxelPool* border_;
    VoxelPool* periodic_;
};

} // ecell4

#endif
