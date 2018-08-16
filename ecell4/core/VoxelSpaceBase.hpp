#ifndef ECELL4_VOXELSPACEBASE_HPP
#define ECELL4_VOXELSPACEBASE_HPP

#include <vector>
#include <set>
#include <map>
#include <stdexcept>
#include <boost/optional.hpp>

#include "Shape.hpp"
// #include "Space.hpp"
#include "Integer3.hpp"
#include "get_mapper_mf.hpp"
#include "Context.hpp"

#ifdef WITH_HDF5
#include "LatticeSpaceHDF5Writer.hpp"
#endif

#include "VoxelPool.hpp"
#include "VacantType.hpp"
#include "ParticleVoxel.hpp"
#include "Particle.hpp"

namespace ecell4
{

#ifdef WIN32_MSC
double rint(const double x);
double round(const double x);
#endif

template <typename T>
static inline std::string get_location_serial(T vp)
{
    if (vp == NULL || vp->location() == NULL || vp->location()->is_vacant()) {
        return "";
    }

    return vp->location()->species().serial();
}


class VoxelSpaceBase
    // : public Space
{
public:

    typedef ParticleVoxel::coordinate_type coordinate_type;

protected:

    typedef utils::get_mapper_mf<Species, VoxelPool*>::type
            voxel_pool_map_type;

    typedef utils::get_mapper_mf<Species, MoleculePool*>::type
            molecule_pool_map_type;

public:

    /*
     * Constructor and Destructor
     */
    VoxelSpaceBase(const Real& voxel_radius) :
        t_(0.0), voxel_radius_(voxel_radius), vacant_(VacantType())
    {}

    VoxelSpaceBase(const Real& voxel_radius, const Shape::dimension_kind& dimension) :
        t_(0.0), voxel_radius_(voxel_radius), vacant_(VacantType(dimension))
    {}

    virtual ~VoxelSpaceBase()
    {
        for (voxel_pool_map_type::iterator itr(voxel_pools_.begin());
             itr != voxel_pools_.end(); ++itr)
        {
            delete itr->second;
        }

        for (molecule_pool_map_type::iterator itr(molecule_pools_.begin());
             itr != molecule_pools_.end(); ++itr)
        {
            delete itr->second;
        }
    }

    /*
     * Static Members
     */

    static inline
    Real
    calculate_voxel_volume(const Real r)
    {
        return 4.0 * sqrt(2.0) * r * r * r;
    }

    static inline
    Real3
    calculate_hcp_lengths(const Real voxel_radius)
    {
        return Real3(
            voxel_radius / sqrt(3.0), // HCP_L
            voxel_radius * sqrt(8.0 / 3.0), // HCP_X
            voxel_radius * sqrt(3.0)); // HCP_Y
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


    /*
     * Space Traits
     */

    const Real t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    /**
     * get the axes lengths of a cuboidal region.
     * @return edge lengths Real3
     */
    virtual const Real3& edge_lengths() const
    {
        throw NotSupported(
            "edge_lengths() is not supported by this space class");
    }

    /**
     * get volume.
     * @return a volume (m^3) Real
     */
    const Real volume() const
    {
        return actual_size() * voxel_volume();
    }

    virtual void save(const std::string& filename) const
    {
        throw NotSupported(
            "save(const std::string) is not supported by this space class");
    }

#ifdef WITH_HDF5
    virtual void save_hdf5(H5::Group* root) const
    {
        throw NotSupported(
            "save_hdf5(H5::Group* root) is not supported by this space class");
    }

    virtual void load_hdf5(const H5::Group& root)
    {
        throw NotSupported(
            "load_hdf5(const H5::Group& root) is not supported by this space class");
    }
#endif


    /*
     * CompartmentSpace Traits
     */
    std::vector<Species> list_species() const;

    bool has_species(const Species& sp) const
    {
        return (voxel_pools_.find(sp) != voxel_pools_.end()
                || molecule_pools_.find(sp) != molecule_pools_.end());
    }

    /*
     * VoxelSpace Traits
     */
    Real voxel_radius() const
    {
        return voxel_radius_;
    }

    Real voxel_volume() const
    {
        return calculate_voxel_volume(voxel_radius_);
    }

    Real get_volume(const Species& sp) const
    {
        return voxel_volume() * num_voxels_exact(sp);
    }

    Real unit_area() const
    {
        const Real r(voxel_radius_);
        return 2.0 * sqrt(3.0) * r * r;
    }

    VoxelPool* vacant() {
        return &vacant_;
    }

    const VoxelPool* vacant() const {
        return &vacant_;
    }

    bool has_voxel(const ParticleID& pid) const;
    Integer num_voxels_exact(const Species& sp) const;
    Integer num_voxels(const Species& sp) const;
    Integer num_voxels() const;

    virtual std::vector<std::pair<ParticleID, ParticleVoxel> > list_voxels() const;
    virtual std::vector<std::pair<ParticleID, ParticleVoxel> > list_voxels(const Species& sp) const;
    virtual std::vector<std::pair<ParticleID, ParticleVoxel> > list_voxels_exact(const Species& sp) const;

    boost::optional<ParticleVoxel> find_voxel(const ParticleID& pid) const;

    std::pair<ParticleID, ParticleVoxel> get_voxel_at(const coordinate_type& coord) const
    {
        const VoxelPool* vp(get_voxel_pool_at(coord));
        return std::make_pair(vp->get_particle_id(coord),
                              ParticleVoxel(
                                  vp->species(),
                                  coord,
                                  vp->radius(),
                                  vp->D(),
                                  get_location_serial(vp)));
    }

    VoxelPool* find_voxel_pool(const Species& sp);
    const VoxelPool* find_voxel_pool(const Species& sp) const;

    bool has_molecule_pool(const Species& sp) const;

    MoleculePool* find_molecule_pool(const Species& sp);
    const MoleculePool* find_molecule_pool(const Species& sp) const;

    virtual VoxelPool* get_voxel_pool_at(const coordinate_type& coord) = 0;
    virtual const VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const = 0;

    /*
     * Coordinate Transformation
     */
    virtual Real3 coordinate2position(const coordinate_type& coord) const = 0;
    virtual coordinate_type position2coordinate(const Real3& pos) const = 0;

    /*
     * Neighbor
     */
    virtual Integer num_neighbors(const coordinate_type& coord) const = 0;

    virtual coordinate_type
    get_neighbor(const coordinate_type& coord, const Integer& nrand) const = 0;

    /*
     * ParticleVoxel Manipulation
     */
    virtual bool update_voxel(const ParticleID& pid, ParticleVoxel v) = 0;
    virtual bool add_voxel(const Species& species, const ParticleID& pid, const coordinate_type& coord) = 0;
    virtual bool remove_voxel(const ParticleID& pid) = 0;
    virtual bool remove_voxel(const coordinate_type& coord) = 0;

    virtual bool can_move(const coordinate_type& src, const coordinate_type& dest) const = 0;

    virtual
    bool
    move(const coordinate_type& src, const coordinate_type& dest,
         const std::size_t candidate=0)
    = 0;

    virtual Integer size() const = 0;
    virtual Integer actual_size() const = 0;

    bool
    make_molecular_type(const Species& sp, Real radius, Real D, const std::string loc);

    bool
    make_structure_type(const Species& sp, Shape::dimension_kind dimension, const std::string loc);

protected:

    VoxelPool* get_voxel_pool(ParticleVoxel v);

    void push_voxels(std::vector<std::pair<ParticleID, ParticleVoxel> >& voxels,
                     const MoleculePool* voxel_pool,
                     const Species& species) const;

protected:

    Real t_;
    Real voxel_radius_;

    VacantType vacant_;

    voxel_pool_map_type    voxel_pools_;
    molecule_pool_map_type molecule_pools_;

};

} // ecell4

#endif /* ECELL4_VOXELSPACEBASE_HPP */
