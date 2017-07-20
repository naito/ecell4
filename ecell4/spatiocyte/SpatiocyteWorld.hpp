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

#include "OneToManyMap.hpp"

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

protected:

    template <typename T>
    class SpaceItem
    {
    public:

        boost::shared_ptr<T> space;
        coordinate_type      offset;
        coordinate_type      inner_offset;

        SpaceItem(T* space, const coordinate_type& offset, const coordinate_type& inner_offset)
            : space(space), offset(offset), inner_offset(inner_offset)
        {}

        SpaceItem(boost::shared_ptr<T>   space,
                  const coordinate_type& offset,
                  const coordinate_type& inner_offset)
            : space(space), offset(offset), inner_offset(inner_offset)
        {}

        operator SpaceItem<const T>() const
        {
            return SpaceItem<const T>(space.get(), offset, inner_offset);
        }
    };

    typedef SpaceItem<VoxelSpaceBase> space_type;
    typedef SpaceItem<const VoxelSpaceBase> const_space_type;

public:

    SpatiocyteWorld(const Real3& edge_lengths,
                    const Real& voxel_radius,
                    const boost::shared_ptr<RandomNumberGenerator>& rng)
        : root_(new default_space_type(edge_lengths, voxel_radius)),
          size_(root_->size()), inner_size_(root_->inner_size()),
          rng_(rng) {}

    SpatiocyteWorld(const Real3& edge_lengths, const Real& voxel_radius) : size_(0), inner_size_(0)
    {
        add_space(new default_space_type(edge_lengths, voxel_radius));
        rng_ = boost::shared_ptr<RandomNumberGenerator>(new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    SpatiocyteWorld(const Real3& edge_lengths = Real3(1, 1, 1))
        : root_(new default_space_type(edge_lengths, edge_lengths[0] / 100.0)),
          size_(root_->size()), inner_size_(root_->inner_size())
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(new GSLRandomNumberGenerator());
        (*rng_).seed();
    }

    SpatiocyteWorld(const std::string filename)
        : root_(new default_space_type(Real3(1, 1, 1), 0.01)),
          size_(root_->size()), inner_size_(root_->inner_size())
    {
        rng_ = boost::shared_ptr<RandomNumberGenerator>(new GSLRandomNumberGenerator());
        this->load(filename);
    }

    SpatiocyteWorld(VoxelSpaceBase* space, const boost::shared_ptr<RandomNumberGenerator>& rng)
        : root_(space), size_(root_->size()), inner_size_(root_->inner_size()), rng_(rng) {}

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

    coordinate_type get_adjoining_or_self(const coordinate_type& coordinate)
    {
        std::vector<coordinate_type> adjoinings(interfaces_.get(coordinate));

        if (adjoinings.empty())
            return coordinate;

        return adjoinings.at(rng_->uniform_int(0, adjoinings.size()-1));
    }


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

    space_type get_corresponding_space(const coordinate_type& coordinate)
    {
        if (root_->size() > coordinate)
            return space_type(root_, 0, 0);

        // should use a binary search algorithm
        for (std::vector<space_type>::iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).offset + (*itr).space->size() > coordinate)
                return *itr;
        }
        throw NotSupported("Out of range");
    }

    const_space_type get_corresponding_space(const coordinate_type& coordinate) const
    {
        if (root_->size() > coordinate)
            return const_space_type(root_, 0, 0);

        // should use a binary search algorithm
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).offset + (*itr).space->size() > coordinate)
                return *itr;
        }
        throw NotSupported("Out of range");
    }

    const_space_type get_corresponding_space_from_inner(const coordinate_type& coordinate) const
    {
        if (root_->inner_size() > coordinate)
            return const_space_type(root_, 0, 0);

        // should use a binary search algorithm
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).inner_offset + (*itr).space->inner_size() > coordinate)
                return *itr;
        }
        throw NotSupported("Out of range");
    }

    space_type get_assigned_space(const ParticleID& pid)
    {
        if (root_->has_particle(pid))
            return space_type(root_, 0, 0);

        for (std::vector<space_type>::iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_particle(pid))
                return *itr;
        }
        throw NotFound("There is no particle having the given ParticleID.");
    }

    const_space_type get_assigned_space(const ParticleID& pid) const
    {
        if (root_->has_particle(pid))
            return const_space_type(root_, 0, 0);

        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_particle(pid))
                return *itr;
        }
        throw NotFound("There is no particle having the given ParticleID.");
    }

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
        return root_->t();
    }

    void set_t(const Real& t)
    {
        root_->set_t(t);
    }

    void save(const std::string& filename) const;
    void load(const std::string& filename);

    /*
     * CompartmentSpaceTraits
     */
    const Real volume() const
    {
        return root_->volume();
    }

    Integer num_species() const
    {
        Integer num_species(root_->num_species());

        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_species += (*itr).space->num_species();
        }

        return num_species;
    }

    bool has_species(const Species& species) const
    {
        if (root_->has_species(species))
            return true;

        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_species(species))
                return true;
        }
        return false;
    }

    // bool has_species_exact(const Species &sp) const;

    Integer num_molecules(const Species& species) const
    {
        Integer num_molecules(root_->num_molecules(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_molecules += (*itr).space->num_molecules(species);
        }
        return num_molecules;
    }

    Integer num_molecules_exact(const Species& species) const
    {
        Integer num_molecules(root_->num_molecules_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_molecules += (*itr).space->num_molecules_exact(species);
        }
        return num_molecules;
    }

    Real get_value(const Species& species) const
    {
        Real value(root_->get_value(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            value += (*itr).space->get_value(species);
        }
        return value;
    }

    Real get_value_exact(const Species& species) const
    {
        Real value(root_->get_value_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
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
        return root_->edge_lengths();
    }

    Real3 actual_lengths() const
    {
        return root_->actual_lengths();
    }

    Integer num_particles() const
    {
        Integer num_particles(root_->num_particles());
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_particles += (*itr).space->num_particles();
        }
        return num_particles;
    }

    Integer num_particles(const Species& species) const
    {
        Integer num_particles(root_->num_particles(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_particles += (*itr).space->num_particles(species);
        }
        return num_particles;
    }

    Integer num_particles_exact(const Species& species) const
    {
        Integer num_particles(root_->num_particles_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_particles += (*itr).space->num_particles_exact(species);
        }
        return num_particles;
    }

    bool has_particle(const ParticleID& pid) const
    {
        if (root_->has_particle(pid))
            return true;

        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_particle(pid))
                return true;
        }
        return false;
    }

    identified_particle get_particle(const ParticleID& pid) const
    {
        if (root_->has_particle(pid))
            return root_->get_particle(pid);

        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_particle(pid))
                return (*itr).space->get_particle(pid);
        }
        throw "Not Found";
    }

    std::vector<identified_particle> list_particles() const
    {
        std::vector<identified_particle> retval(root_->list_particles());
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_particle> particles((*itr).space->list_particles());
            retval.insert(retval.end(), particles.begin(), particles.end());
        }
        return retval;
    }

    std::vector<identified_particle> list_particles(const Species& species) const
    {
        std::vector<identified_particle> retval(root_->list_particles(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_particle> particles((*itr).space->list_particles(species));
            retval.insert(retval.end(), particles.begin(), particles.end());
        }
        return retval;
    }

    std::vector<identified_particle> list_particles_exact(const Species& species) const
    {
        std::vector<identified_particle> retval(root_->list_particles_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
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
        const_space_type item(get_corresponding_space(coordinate));
        return item.space->coordinate2position(coordinate - item.offset);
    }

    coordinate_type position2coordinate(const Real3& position) const
    {
        return root_->position2coordinate(position);
    }

    Integer num_neighbors(const coordinate_type& coordinate) const
    {
        const_space_type item(get_corresponding_space(coordinate));
        return item.space->num_neighbors(coordinate - item.offset);
    }

    coordinate_type get_neighbor(coordinate_type coordinate, Integer nrand) const
    {
        const_space_type item(get_corresponding_space(coordinate));
        return item.space->get_neighbor(coordinate - item.offset, nrand) + item.offset;
    }

    coordinate_type get_neighbor_boundary(coordinate_type coordinate, Integer nrand) const
    {
        const_space_type item(get_corresponding_space(coordinate));
        return item.space->get_neighbor_boundary(coordinate - item.offset, nrand) + item.offset;
    }

    Integer size() const
    {
        return size_;
    }

    const Integer3 shape() const
    {
        return root_->shape();
    }

    const Integer inner_size() const
    {
        return inner_size_;
    }

    //TODO
    coordinate_type inner2coordinate(const coordinate_type coordinate) const
    {
        const_space_type item(get_corresponding_space_from_inner(coordinate));
        return item.space->inner2coordinate(coordinate - item.inner_offset);
    }

    /*
     * VoxelSpaceTraits
     */
    Real voxel_radius() const
    {
        return root_->voxel_radius();
    }

    Real unit_voxel_volume() const
    {
        return root_->unit_voxel_volume();
    }

    Real voxel_volume() const
    {
        return root_->voxel_volume();
    }

    Real unit_area() const
    {
        return root_->unit_area();
    }

    std::vector<Species> list_species() const
    {
        std::vector<Species> species(root_->list_species());
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<Species> new_species((*itr).space->list_species());
            species.insert(species.end(), new_species.begin(), new_species.end());
        }
        return species;
    }

    Integer num_voxels() const
    {
        Integer num_voxels(root_->num_voxels());
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_voxels += (*itr).space->num_voxels();
        }
        return num_voxels;
    }

    Integer num_voxels(const Species& species) const
    {
        Integer num_voxels(root_->num_voxels(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_voxels += (*itr).space->num_voxels(species);
        }
        return num_voxels;
    }

    Integer num_voxels_exact(const Species& species) const
    {
        Integer num_voxels(root_->num_voxels_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            num_voxels += (*itr).space->num_voxels_exact(species);
        }
        return num_voxels;
    }

    bool has_voxel(const ParticleID& pid) const
    {
        if (root_->has_voxel(pid))
            return true;

        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_voxel(pid))
                return true;
        }
        return false;
    }

    identified_voxel get_voxel(const ParticleID& pid) const
    {
        const_space_type item(get_assigned_space(pid));
        identified_voxel retval(item.space->get_voxel(pid));
        retval.second.coordinate() -= item.offset;
        return retval;
    }

    std::vector<identified_voxel> list_voxels() const
    {
        std::vector<identified_voxel> voxels(root_->list_voxels());
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_voxel> new_voxels((*itr).space->list_voxels());
            voxels.insert(voxels.end(), new_voxels.begin(), new_voxels.end());
        }
        return voxels;
    }

    std::vector<identified_voxel> list_voxels(const Species& species) const
    {
        std::vector<identified_voxel> voxels(root_->list_voxels(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_voxel> new_voxels((*itr).space->list_voxels(species));
            voxels.insert(voxels.end(), new_voxels.begin(), new_voxels.end());
        }
        return voxels;
    }

    std::vector<identified_voxel> list_voxels_exact(const Species& species) const
    {
        std::vector<identified_voxel> voxels(root_->list_voxels_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            std::vector<identified_voxel> new_voxels((*itr).space->list_voxels_exact(species));
            voxels.insert(voxels.end(), new_voxels.begin(), new_voxels.end());
        }
        return voxels;
    }

    // bool has_voxel_pool(const Species& sp) const;

    VoxelPool* find_voxel_pool(const Species& species)
    {
        if (root_->has_species(species))
            return root_->find_voxel_pool(species);

        for (std::vector<space_type>::iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_species(species))
                return (*itr).space->find_voxel_pool(species);
        }
        throw "Not Found";
    }

    const VoxelPool* find_voxel_pool(const Species& species) const
    {
        if (root_->has_species(species))
            return root_->find_voxel_pool(species);

        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_species(species))
                return (*itr).space->find_voxel_pool(species);
        }
        throw "Not Found";
    }

    bool has_molecule_pool(const Species& species) const
    {
        if (root_->has_molecule_pool(species))
            return true;

        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_molecule_pool(species))
                return true;
        }
        return false;
    }

    MoleculePool* find_molecule_pool(const Species& species)
    {
        if (root_->has_molecule_pool(species))
            return root_->find_molecule_pool(species);

        for (std::vector<space_type>::iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_molecule_pool(species))
                return (*itr).space->find_molecule_pool(species);
        }
        throw "Not Found";
    }

    const MoleculePool* find_molecule_pool(const Species& species) const
    {
        if (root_->has_molecule_pool(species))
            return root_->find_molecule_pool(species);

        for (std::vector<space_type>::const_iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).space->has_molecule_pool(species))
                return (*itr).space->find_molecule_pool(species);
        }
        throw "Not Found";
    }

    bool remove_particle(const ParticleID& pid)
    {
        space_type item(get_assigned_space(pid));
        return item.space->remove_particle(pid);
    }

    bool on_structure(const Voxel& voxel)
    {
        space_type item(get_corresponding_space(voxel.coordinate()));
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
        const_space_type item(get_corresponding_space(coordinate));
        return item.space->get_voxel_at(coordinate - item.offset);
    }

    const VoxelPool* get_voxel_pool_at(const coordinate_type& coordinate) const
    {
        const_space_type item(get_corresponding_space(coordinate));
        return item.space->get_voxel_pool_at(coordinate - item.offset);
    }

    /* Not wrapped function
    virtual const Particle particle_at(const coordinate_type& coord) const;
    */

    bool update_voxel(const ParticleID& pid, const Voxel& v)
    {
        // XXX
        return root_->update_voxel(pid, v);
    }

    bool remove_voxel(const ParticleID pid)
    {
        space_type item(get_assigned_space(pid));
        return item.space->remove_voxel(pid);
    }

    bool remove_voxel(const coordinate_type& coordinate)
    {
        space_type item(get_corresponding_space(coordinate));
        return item.space->remove_voxel(coordinate - item.offset);
    }

    bool can_move(const coordinate_type& src, const coordinate_type& dest) const
    {
        // XXX
        return root_->can_move(src, dest);
    }

    bool move(const coordinate_type& src, const coordinate_type& dest,
              const std::size_t candidate=0)
    {
        // XXX
        return root_->move(src, dest, candidate);
    }

    std::pair<coordinate_type, bool>
    move_to_neighbor(VoxelPool* const& from_mt, VoxelPool* const& loc,
                     coordinate_id_pair_type& info, const Integer nrand)
    {
        // XXX
        return root_->move_to_neighbor(from_mt, loc, info, nrand);
    }

protected:

    boost::shared_ptr<VoxelSpaceBase> root_;
    Integer size_;
    Integer inner_size_;
    std::vector<space_type> spaces_;

    OneToManyMap<coordinate_type> interfaces_;
    OneToManyMap<coordinate_type> neighbors_;

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
