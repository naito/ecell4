#ifndef __ECELL4_VOXELSPACEBASE_HPP
#define __ECELL4_VOXELSPACEBASE_HPP

#include <vector>
#include <set>
#include <map>

#include "get_mapper_mf.hpp"
#include "Space.hpp"
#include "Voxel.hpp"
#include "VoxelPool.hpp"
#include "Shape.hpp"
#include "Integer3.hpp"

#ifdef WITH_HDF5
#include "LatticeSpaceHDF5Writer.hpp"
#endif

namespace ecell4
{

class VoxelSpaceBase : public Space
{
public:

    typedef Voxel::coordinate_type             coordinate_type;
    typedef VoxelPool::coordinate_id_pair_type coordinate_id_pair_type;
    typedef std::pair<ParticleID, Voxel>       identified_voxel;

protected:

    typedef utils::get_mapper_mf<Species, boost::shared_ptr<VoxelPool> >::type
            voxel_pool_map_type;

    typedef utils::get_mapper_mf<Species, boost::shared_ptr<MoleculePool> >::type
            molecule_pool_map_type;

public:

    VoxelSpaceBase(const Real& voxel_radius);
    virtual ~VoxelSpaceBase();

    /*
     * for LatticeSpace
     */

    virtual Real3 coordinate2position(const coordinate_type& coord) const = 0;
    virtual coordinate_type position2coordinate(const Real3& pos) const = 0;

    virtual Integer num_neighbors(const coordinate_type& coord) const = 0;
    virtual coordinate_type
        get_neighbor(const coordinate_type& coord, const Integer& nrand) const = 0;
    virtual coordinate_type
        get_neighbor_boundary(const coordinate_type& coord, const Integer& nrand) const = 0;

    virtual Integer size() const = 0;
    virtual Integer3 shape() const = 0;

    // For spaces having the inner structure
    virtual coordinate_type inner2coordinate(const coordinate_type inner) const = 0;
    virtual Integer inner_size() const = 0;

    /*
     * VoxelSpaceTraits
     */
    inline Real voxel_radius() const
    {
        return voxel_radius_;
    }

    inline Real voxel_volume() const
    {
        return unit_voxel_volume() * pow(voxel_radius_, 3);
    }

    inline Real unit_area() const
    {
        return 2.0 * sqrt(3.0) * pow(voxel_radius_, 2);
    }

    std::vector<Species> list_species() const;

    Integer num_voxels() const;
    Integer num_voxels(const Species& sp) const;
    Integer num_voxels_exact(const Species& sp) const;
    bool has_voxel(const ParticleID& pid) const;
    identified_voxel get_voxel(const ParticleID& pid) const;
    std::vector<identified_voxel> list_voxels() const;
    std::vector<identified_voxel> list_voxels(const Species& sp) const;
    std::vector<identified_voxel> list_voxels_exact(const Species& sp) const;

    bool has_voxel_pool(const Species& sp) const;
    const VoxelPool* find_voxel_pool(const Species& sp) const;
          VoxelPool* find_voxel_pool(const Species& sp);

    bool has_molecule_pool(const Species& sp) const;
    const MoleculePool* find_molecule_pool(const Species& sp) const;
          MoleculePool* find_molecule_pool(const Species& sp);

    bool remove_particle(const ParticleID& pid);

    bool on_structure(const Voxel& v);
    bool make_structure_type(const Species& sp,
                             Shape::dimension_kind dimension,
                             const std::string loc);
    bool make_interface_type(const Species& sp,
                             Shape::dimension_kind dimension,
                             const std::string loc);
    bool make_molecular_pool(const Species& sp,
                             Real radius,
                             Real D,
                             const std::string loc);

    virtual Real unit_voxel_volume() const = 0;

    virtual identified_voxel get_voxel_at(const coordinate_type& coord) const = 0;
    virtual const VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const = 0;
    virtual       VoxelPool* get_voxel_pool_at(const coordinate_type& coord)       = 0;
    virtual const Particle particle_at(const coordinate_type& coord) const = 0;

    virtual bool update_voxel(const ParticleID& pid, const Voxel& v) = 0;
    virtual bool remove_voxel(const ParticleID& pid) = 0;
    virtual bool remove_voxel(const coordinate_type& coord) = 0;

    virtual bool can_move(const coordinate_type& src, const coordinate_type& dest) const = 0;
    virtual bool
        move(const coordinate_type& src,
             const coordinate_type& dest,
             const std::size_t candidate = 0)
        = 0;
    virtual std::pair<coordinate_type, bool>
        move_to_neighbor(VoxelPool* const& from,
                         VoxelPool* const& loc,
                         coordinate_id_pair_type& info,
                         const Integer nrand)
        = 0;

protected:

    virtual Integer count_voxels(const boost::shared_ptr<VoxelPool>& vp) const = 0;
    VoxelPool* get_voxel_pool(const Voxel& voxel);

private:
    std::string get_location_serial(const boost::shared_ptr<VoxelPool>& voxel_pool) const;
    void push_voxels(std::vector<identified_voxel>& voxels,
                     const boost::shared_ptr<MoleculePool>& voxel_pool) const;
    void push_voxels(std::vector<identified_voxel>& voxels,
                     const boost::shared_ptr<VoxelPool>& voxel_pool) const;
    VoxelPool* get_vp_from_serial(const std::string& serial);

    template<typename T>
    bool make_location_type(const Species& sp,
                            Shape::dimension_kind dimension,
                            const std::string loc);

public:

    /*
     * SpaceTraits
     */

    const Real t() const;
    void set_t(const Real& t);

    virtual void save(const std::string& filename) const
    {
        throw NotSupported(
            "save(const std::string) is not supported by this space class");
    }

#ifdef WITH_HDF5
    virtual void save_hdf5(H5::Group* root) const
    {
        throw NotSupported(
            "load(H5::Group* root) is not supported by this space class");
    }

    virtual void load_hdf5(const H5::Group& root)
    {
        throw NotSupported(
            "load(const H5::Group& root) is not supported by this space class");
    }
#endif

    /*
     * CompartmentSpaceTraits
     */
    const Real volume() const;
    Integer num_species() const;
    bool has_species(const Species& sp) const;

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;

    Real get_value(const Species& sp) const;
    Real get_value_exact(const Species& sp) const;

    /*
     * ParticleSpaceTraits
     */
    // virtual Real3 edge_lengths() const = 0;
    // virtual Real3 actual_lengths() const = 0;

    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;

    bool has_particle(const ParticleID& pid) const;

    identified_particle get_particle(const ParticleID& pid) const;
    std::vector<identified_particle> list_particles() const;
    std::vector<identified_particle> list_particles(const Species& sp) const;
    std::vector<identified_particle> list_particles_exact(const Species& sp) const;


protected:

    Real t_;
    Real voxel_radius_;
    VoxelPool* vacant_;

    voxel_pool_map_type voxel_pools_;
    molecule_pool_map_type molecule_pools_;

};

} // ecell4

#endif /* __ECELL4_VOXELSPACEBASE_HPP */
