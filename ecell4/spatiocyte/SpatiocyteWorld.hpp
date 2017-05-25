#ifndef __ECELL4_LATTICE_LATTICE_WORLD_HPP
#define __ECELL4_LATTICE_LATTICE_WORLD_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <ecell4/core/VoxelPool.hpp>
#include <ecell4/core/VoxelSpaceBase.hpp>
#include <ecell4/core/LatticeSpaceVectorImpl.hpp>
#include <ecell4/core/OffLatticeSpace.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Shape.hpp>

namespace ecell4
{

namespace spatiocyte
{

struct MoleculeInfo
{
    const Real radius;
    const Real D;
    const std::string loc;
};

class SpatiocyteWorld : public Space
{
public:

    typedef LatticeSpaceVectorImpl default_space_type;

    typedef VoxelSpaceBase::coordinate_id_pair_type coordinate_id_pair_type;
    typedef VoxelSpaceBase::coordinate_type coordinate_type;

    typedef VoxelSpaceBase::identified_voxel identified_voxel;

    struct SpaceItem
    {
        boost::shared_ptr<VoxelSpaceBase> space;
        coordinate_type                   offset;
        coordinate_type                   inner_offset;
    };

public:

    SpatiocyteWorld(const Real3& edge_lengths, const Real& voxel_radius,
        const boost::shared_ptr<RandomNumberGenerator>& rng)
        : size_(0), inner_size_(0), rng_(rng)
    {
        add_space(new default_space_type(edge_lengths, voxel_radius));
    }

    SpatiocyteWorld(const Real3& edge_lengths, const Real& voxel_radius) : size_(0), inner_size_(0)
    {
        add_space(new default_space_type(edge_lengths, voxel_radius));
        rng_ = boost::shared_ptr<RandomNumberGenerator>(new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    SpatiocyteWorld(const Real3& edge_lengths = Real3(1, 1, 1)) : size_(0), inner_size_(0)
    {
        add_space(new default_space_type(edge_lengths, edge_lengths[0] / 100.0));
        rng_ = boost::shared_ptr<RandomNumberGenerator>(new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    SpatiocyteWorld(const std::string filename) : size_(0), inner_size_(0)
    {
        add_space(new default_space_type(Real3(1, 1, 1), 0.01));
        rng_ = boost::shared_ptr<RandomNumberGenerator>(new GSLRandomNumberGenerator());
        this->load(filename);
    }

    SpatiocyteWorld(VoxelSpaceBase* space, const boost::shared_ptr<RandomNumberGenerator>& rng)
        : size_(0), inner_size_(0), rng_(rng)
    {
        add_space(space);
    }

    /*
     * Class functions
     */

    /**
     * draw attributes of species and return it as a molecule info.
     * @param sp a species
     * @return info a molecule info
     */
    MoleculeInfo get_molecule_info(const Species& sp) const;

    std::pair<identified_voxel, bool> new_voxel(const Voxel& v)
    {
        ParticleID pid(sidgen_());
        const bool is_succeeded(update_voxel(pid, v));
        return std::make_pair(std::make_pair(pid, v), is_succeeded);
    }

    std::pair<identified_voxel, bool> new_voxel(const Species& sp, const coordinate_type& coord)
    {
        const MoleculeInfo minfo(get_molecule_info(sp));
        return new_voxel(Voxel(sp, coord, minfo.radius, minfo.D, minfo.loc));
    }

    std::pair<coordinate_type, bool>
    check_neighbor(const coordinate_type coord, const std::string& loc);

    boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    identified_voxel
    make_pid_voxel_pair(const VoxelPool* vpool, const coordinate_id_pair_type& info) const
    {
        const std::string loc(vpool->location()->is_vacant() ?
                "" : vpool->location()->species().serial());
        return identified_voxel(
            ParticleID(info.pid),
            Voxel(vpool->species(), info.coordinate, vpool->radius(), vpool->D(), loc));
    }

    identified_voxel
    make_pid_voxel_pair(const VoxelPool* vpool, const coordinate_type& coord) const
    {
        const ParticleID pid(vpool->get_particle_id(coord));
        const coordinate_id_pair_type info(pid, coord);
        return make_pid_voxel_pair(vpool, info);
    }

    identified_voxel choice(const Species& sp)
    {
        const MoleculePool* mpool(find_molecule_pool(sp));
        const Integer i(rng_->uniform_int(0, mpool->size()-1));
        return make_pid_voxel_pair(mpool, mpool->at(i));
    }

    boost::shared_ptr<Model> lock_model() const
    {
        return model_.lock();
    }

    // XXX: Not implemented
    const Species& draw_species(const Species& pttrn) const;

    // XXX: Never called
    Integer add_neighbors(const Species& sp, const coordinate_type center);


public:

    /*
     * Python API
     */

    void set_value(const Species& sp, const Real value);

    /**
     * create and add a new particle
     * @param p a particle
     * @return a pair of a pair of pid (a particle id) and p (a particle)
     * and bool (if it's succeeded or not)
     */
    std::pair<identified_particle, bool> new_particle(const Particle& p);
    std::pair<identified_particle, bool> new_particle(const Species& sp, const Real3& pos)
    {
        const MoleculeInfo info(get_molecule_info(sp));
        return new_particle(Particle(sp, pos, info.radius, info.D));
    }

    bool update_particle(const ParticleID& pid, const Particle& p);

    std::vector<identified_particle> list_structure_particles() const;
    std::vector<identified_particle> list_non_structure_particles() const;

    std::pair<identified_voxel, bool>
    new_voxel_structure(const Species& sp, const coordinate_type& coord);

    std::pair<identified_voxel, bool>
    new_voxel_interface(const Species& sp, const coordinate_type& coord);

    Integer add_structure(const Species& sp, const boost::shared_ptr<const Shape> shape);
    Integer add_interface(const Species& sp);

    std::vector<Species> list_structure_species() const;
    std::vector<Species> list_non_structure_species() const;

    bool add_molecules(const Species& sp, const Integer& num);
    bool add_molecules(const Species& sp, const Integer& num,
                       const boost::shared_ptr<const Shape> shape);
    void remove_molecules(const Species& sp, const Integer& num);

    void bind_to(boost::shared_ptr<Model> model);


protected:

    void add_space(VoxelSpaceBase *space);
    SpaceItem& get_corresponding_space(const coordinate_type& coordinate);
    const SpaceItem& get_corresponding_space(const coordinate_type& coordinate) const;
    SpaceItem& get_corresponding_space_from_inner(const coordinate_type& coordinate);
    const SpaceItem& get_corresponding_space_from_inner(const coordinate_type& coordinate) const;
    SpaceItem& get_assigned_space(const ParticleID& pid);
    const SpaceItem& get_assigned_space(const ParticleID& pid) const;

    std::pair<identified_voxel, bool> new_voxel_structure(const Voxel& v);
    Integer add_structure3(const Species& sp, const boost::shared_ptr<const Shape> shape);
    Integer add_structure2(const Species& sp, const boost::shared_ptr<const Shape> shape);
    bool is_surface_voxel(const coordinate_type coord,
            const boost::shared_ptr<const Shape> shape) const;

public:

    /* Wrapped functions */

    /*
     * SpaceTraits
     */
    const Real t() const
    {
        if (spaces_.size() == 0)
            return 0.0;

        return spaces_.at(0).space->t();
    }

    void set_t(const Real& t)
    {
        if (spaces_.size() != 0)
            spaces_.at(0).space->set_t(t);
    }

    void save(const std::string& filename) const;
    void load(const std::string& filename);

    /*
     * CompartmentSpaceTraits
     */
    const Real volume() const
    {
        if (spaces_.size() == 0)
            return 0.0;

        return spaces_.at(0).space->volume();
    }

    Integer num_species() const
    {
        Integer num_species = 0;

        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_species += (*itr).space->num_species();
        }

        return num_species;
    }

    bool has_species(const Species& species) const
    {
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_species(species))
                return true;
        }
        return false;
    }

    // bool has_species_exact(const Species &sp) const;

    Integer num_molecules(const Species& species) const
    {
        Integer num_molecules = 0;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_molecules += (*itr).space->num_molecules(species);
        }
        return num_molecules;
    }

    Integer num_molecules_exact(const Species& species) const
    {
        Integer num_molecules = 0;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_molecules += (*itr).space->num_molecules_exact(species);
        }
        return num_molecules;
    }

    Real get_value(const Species& species) const
    {
        Real value = 0;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            value += (*itr).space->get_value(species);
        }
        return value;
    }

    Real get_value_exact(const Species& species) const
    {
        Real value = 0;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            value += (*itr).space->get_value_exact(species);
        }
        return value;
    }


    /*
     * ParticleSpaceTraits
     */
    const Real3& edge_lengths() const
    {
        // This code returns the reference to a stack memory
        // if (spaces_.size() == 0)
        //     return Real3();

        return spaces_.at(0).space->edge_lengths();
    }

    Real3 actual_lengths() const
    {
        if (spaces_.size() == 0)
            return Real3();

        return spaces_.at(0).space->actual_lengths();
    }

    Integer num_particles() const
    {
        Integer num_particles(0);
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_particles += (*itr).space->num_particles();
        }
        return num_particles;
    }

    Integer num_particles(const Species& species) const
    {
        Integer num_particles(0);
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_particles += (*itr).space->num_particles(species);
        }
        return num_particles;
    }

    Integer num_particles_exact(const Species& species) const
    {
        Integer num_particles(0);
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_particles += (*itr).space->num_particles_exact(species);
        }
        return num_particles;
    }

    bool has_particle(const ParticleID& pid) const
    {
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_particle(pid))
                return true;
        }
        return false;
    }

    identified_particle get_particle(const ParticleID& pid) const
    {
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_particle(pid))
                return (*itr).space->get_particle(pid);
        }
        throw "Not Found";
    }

    std::vector<identified_particle> list_particles() const
    {
        std::vector<identified_particle> retval;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_particle> particles((*itr).space->list_particles());
            retval.insert(retval.end(), particles.begin(), particles.end());
        }
        return retval;
    }

    std::vector<identified_particle> list_particles(const Species& species) const
    {
        std::vector<identified_particle> retval;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_particle> particles((*itr).space->list_particles(species));
            retval.insert(retval.end(), particles.begin(), particles.end());
        }
        return retval;
    }

    std::vector<identified_particle> list_particles_exact(const Species& species) const
    {
        std::vector<identified_particle> retval;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_particle> particles((*itr).space->list_particles_exact(species));
            retval.insert(retval.end(), particles.begin(), particles.end());
        }
        return retval;
    }

    /*
     * LatticeSpaceTraits
     */
    const Real3 coordinate2position(const coordinate_type& coordinate) const
    {
        const SpaceItem& item(get_corresponding_space(coordinate));
        return item.space->coordinate2position(coordinate - item.offset);
    }

    coordinate_type position2coordinate(const Real3& position) const
    {
        return spaces_.at(0).space->position2coordinate(position);
    }

    Integer num_neighbors(const coordinate_type& coordinate) const
    {
        const SpaceItem& item(get_corresponding_space(coordinate));
        return item.space->num_neighbors(coordinate - item.offset);
    }

    coordinate_type get_neighbor(coordinate_type coordinate, Integer nrand) const
    {
        const SpaceItem& item(get_corresponding_space(coordinate));
        return item.space->get_neighbor(coordinate - item.offset, nrand) + item.offset;
    }

    coordinate_type get_neighbor_boundary(coordinate_type coordinate, Integer nrand) const
    {
        const SpaceItem& item(get_corresponding_space(coordinate));
        return item.space->get_neighbor_boundary(coordinate - item.offset, nrand) + item.offset;
    }

    Integer size() const
    {
        return size_;
    }

    const Integer3 shape() const
    {
        if (spaces_.size() == 0)
            return Integer3();

        return spaces_.at(0).space->shape();
    }

    const Integer inner_size() const
    {
        return inner_size_;
    }

    //TODO
    coordinate_type inner2coordinate(const coordinate_type coordinate) const
    {
        const SpaceItem& item(get_corresponding_space_from_inner(coordinate));
        return item.space->inner2coordinate(coordinate - item.inner_offset);
    }

    /*
     * VoxelSpaceTraits
     */
    Real voxel_radius() const
    {
        if (spaces_.size() == 0)
            return 0.0;

        return spaces_.at(0).space->voxel_radius();
    }

    Real unit_voxel_volume() const
    {
        if (spaces_.size() == 0)
            return 0.0;

        return spaces_.at(0).space->unit_voxel_volume();
    }

    Real voxel_volume() const
    {
        if (spaces_.size() == 0)
            return 0.0;

        return spaces_.at(0).space->voxel_volume();
    }

    Real unit_area() const
    {
        if (spaces_.size() == 0)
            return 0.0;

        return spaces_.at(0).space->unit_area();
    }

    std::vector<Species> list_species() const
    {
        std::vector<Species> species;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<Species> new_species((*itr).space->list_species());
            species.insert(species.end(), new_species.begin(), new_species.end());
        }
        return species;
    }

    Integer num_voxels() const
    {
        Integer num_voxels = 0;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_voxels += (*itr).space->num_voxels();
        }
        return num_voxels;
    }

    Integer num_voxels(const Species& species) const
    {
        Integer num_voxels = 0;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_voxels += (*itr).space->num_voxels(species);
        }
        return num_voxels;
    }

    Integer num_voxels_exact(const Species& species) const
    {
        Integer num_voxels = 0;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_voxels += (*itr).space->num_voxels_exact(species);
        }
        return num_voxels;
    }

    bool has_voxel(const ParticleID& pid) const
    {
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_voxel(pid))
                return true;
        }
        return false;
    }

    identified_voxel get_voxel(const ParticleID& pid) const
    {
        const SpaceItem& item(get_assigned_space(pid));
        identified_voxel retval(item.space->get_voxel(pid));
        retval.second.coordinate() -= item.offset;
        return retval;
    }

    std::vector<identified_voxel> list_voxels() const
    {
        std::vector<identified_voxel> voxels;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_voxel> new_voxels((*itr).space->list_voxels());
            voxels.insert(voxels.end(), new_voxels.begin(), new_voxels.end());
        }
        return voxels;
    }

    std::vector<identified_voxel> list_voxels(const Species& species) const
    {
        std::vector<identified_voxel> voxels;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_voxel> new_voxels((*itr).space->list_voxels(species));
            voxels.insert(voxels.end(), new_voxels.begin(), new_voxels.end());
        }
        return voxels;
    }

    std::vector<identified_voxel> list_voxels_exact(const Species& species) const
    {
        std::vector<identified_voxel> voxels;
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_voxel> new_voxels((*itr).space->list_voxels_exact(species));
            voxels.insert(voxels.end(), new_voxels.begin(), new_voxels.end());
        }
        return voxels;
    }

    // bool has_voxel_pool(const Species& sp) const;

    VoxelPool* find_voxel_pool(const Species& species)
    {
        for (std::vector<SpaceItem>::iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_species(species))
                return (*itr).space->find_voxel_pool(species);
        }
        throw "Not Found";
    }

    const VoxelPool* find_voxel_pool(const Species& species) const
    {
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_species(species))
                return (*itr).space->find_voxel_pool(species);
        }
        throw "Not Found";
    }

    bool has_molecule_pool(const Species& species) const
    {
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_molecule_pool(species))
                return true;
        }
        return false;
    }

    MoleculePool* find_molecule_pool(const Species& species)
    {
        for (std::vector<SpaceItem>::iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_molecule_pool(species))
                return (*itr).space->find_molecule_pool(species);
        }
        throw "Not Found";
    }

    const MoleculePool* find_molecule_pool(const Species& species) const
    {
        for (std::vector<SpaceItem>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_molecule_pool(species))
                return (*itr).space->find_molecule_pool(species);
        }
        throw "Not Found";
    }

    bool remove_particle(const ParticleID& pid)
    {
        SpaceItem& item(get_assigned_space(pid));
        return item.space->remove_particle(pid);
    }

    bool on_structure(const Voxel& voxel)
    {
        SpaceItem& item(get_corresponding_space(voxel.coordinate()));
        const Voxel modified(voxel.species(),
                             voxel.coordinate() - item.offset,
                             voxel.radius(),
                             voxel.D(),
                             voxel.loc());
        return item.space->on_structure(modified);
    }

    /* Not wrapped functions
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
     */

    identified_voxel get_voxel_at(const coordinate_type& coordinate) const
    {
        const SpaceItem& item(get_corresponding_space(coordinate));
        return item.space->get_voxel_at(coordinate - item.offset);
    }

    const VoxelPool* get_voxel_pool_at(const coordinate_type& coordinate) const
    {
        const SpaceItem& item(get_corresponding_space(coordinate));
        return item.space->get_voxel_pool_at(coordinate - item.offset);
    }

    /* Not wrapped function
    virtual const Particle particle_at(const coordinate_type& coord) const;
    */

    bool update_voxel(const ParticleID& pid, const Voxel& v)
    {
        return spaces_.at(0).space->update_voxel(pid, v);
    }

    bool remove_voxel(const ParticleID pid)
    {
        SpaceItem& item(get_assigned_space(pid));
        return item.space->remove_voxel(pid);
    }

    bool remove_voxel(const coordinate_type& coordinate)
    {
        SpaceItem& item(get_corresponding_space(coordinate));
        return item.space->remove_voxel(coordinate - item.offset);
    }

    bool can_move(const coordinate_type& src, const coordinate_type& dest) const
    {
        return (*spaces_.at(0).space).can_move(src, dest);
    }

    bool move(const coordinate_type& src, const coordinate_type& dest,
              const std::size_t candidate=0)
    {
        return spaces_.at(0).space->move(src, dest, candidate);
    }

    std::pair<coordinate_type, bool>
    move_to_neighbor(VoxelPool* const& from_mt, VoxelPool* const& loc,
                     coordinate_id_pair_type& info, const Integer nrand)
    {
        return spaces_.at(0).space->move_to_neighbor(from_mt, loc, info, nrand);
    }

protected:

    Integer size_;
    Integer inner_size_;
    std::vector<SpaceItem> spaces_;

    boost::shared_ptr<RandomNumberGenerator> rng_;
    SerialIDGenerator<ParticleID> sidgen_;

    boost::weak_ptr<Model> model_;
};

SpatiocyteWorld* create_spatiocyte_world_cell_list_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const Integer3& matrix_sizes,
    const boost::shared_ptr<RandomNumberGenerator>& rng);
SpatiocyteWorld* create_spatiocyte_world_vector_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng);
SpatiocyteWorld* create_spatiocyte_world_offlattice_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng);

/**
 * Alias functions for Cython
 */

inline SpatiocyteWorld* create_spatiocyte_world_cell_list_impl_alias(
    const Real3& edge_lengths, const Real& voxel_radius,
    const Integer3& matrix_sizes,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return create_spatiocyte_world_cell_list_impl(
        edge_lengths, voxel_radius, matrix_sizes, rng);
}

inline SpatiocyteWorld* create_spatiocyte_world_vector_impl_alias(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return create_spatiocyte_world_vector_impl(edge_lengths, voxel_radius, rng);
}

inline SpatiocyteWorld* create_spatiocyte_world_offlattice_impl_alias(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return create_spatiocyte_world_offlattice_impl(edge_lengths, voxel_radius, rng);
}

} // spatiocyte

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_WORLD_HPP */
