#ifndef ECELL4_LATTICE_LATTICE_WORLD_HPP
#define ECELL4_LATTICE_LATTICE_WORLD_HPP

#include <sstream>
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

#include "OffsetSpace.hpp"
#include "OneToManyMap.hpp"
#include "utils.hpp"

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

    typedef OffsetSpace<VoxelSpaceBase> space_type;
    typedef OffsetSpace<const VoxelSpaceBase> const_space_type;

public:

    SpatiocyteWorld(const Real3& edge_lengths,
                    const Real& voxel_radius,
                    const boost::shared_ptr<RandomNumberGenerator>& rng)
        : root_(new default_space_type(edge_lengths, voxel_radius)),
          size_(root_->size()), inner_size_(root_->inner_size()),
          rng_(rng) {}

    SpatiocyteWorld(const Real3& edge_lengths, const Real& voxel_radius)
        : root_(new default_space_type(edge_lengths, voxel_radius)),
          size_(root_->size()), inner_size_(root_->inner_size())
    {
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

    boost::optional<coordinate_type>
    check_neighbor(const coordinate_type coord, const std::string& loc);

    coordinate_type
    pick_neighbor(const coordinate_type& coordinate)
    {
        const_space_type item(get_space(coordinate));
        const std::size_t num_neighbors(item.num_neighbors(coordinate));
        assert(num_neighbors != 0);
        if (num_neighbors == 1)
            return item.get_neighbor_boundary(coordinate, 0);
        return item.get_neighbor_boundary(coordinate, rng()->uniform_int(1, num_neighbors-1));
    }

    boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

    identified_voxel
    make_pid_voxel_pair(const VoxelPool* vpool, const coordinate_id_pair_type& info) const
    {
        return identified_voxel(
            ParticleID(info.pid),
            Voxel(vpool->species(),
                  info.coordinate,
                  vpool->radius(),
                  vpool->D(),
                  vpool->get_location_serial()));
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
        std::pair<const MoleculePool*, coordinate_type> mpool(find_molecule_pool(sp));

        if (mpool.first->size() == 0)
        {
            throw NotFound("There is no corresponding voxel.");
        }
        const Integer i(mpool.first->size() == 1 ? 0 : rng_->uniform_int(0, mpool.first->size()-1));

        coordinate_id_pair_type info(mpool.first->at(i));
        info.coordinate += mpool.second;
        return make_pid_voxel_pair(mpool.first, info);
    }

    boost::shared_ptr<Model> lock_model() const
    {
        return model_.lock();
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

    void add_space(VoxelSpaceBase *space);

protected:

    space_type get_space_mut(const coordinate_type& coordinate)
    {
        if (root_->size() > coordinate)
            return space_type(root_, 0, 0);

        // should use a binary search algorithm
        for (std::vector<space_type>::iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).has(coordinate))
                return *itr;
        }
        throw NotSupported("Out of range");
    }

    const_space_type get_space(const coordinate_type& coordinate) const
    {
        if (root_->size() > coordinate)
            return const_space_type(root_, 0, 0);

        // should use a binary search algorithm
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has(coordinate))
                return *itr;
        }
        throw NotSupported("Out of range");
    }

    const_space_type get_space_from_inner(const coordinate_type& coordinate) const
    {
        if (root_->inner_size() > coordinate)
            return const_space_type(root_, 0, 0);

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
            if ((*itr).has_inner(coordinate))
                return *itr;

        throw NotSupported("Out of range");
    }

    space_type get_space_mut(const ParticleID& pid)
    {
        if (root_->has_particle(pid))
            return space_type(root_, 0, 0);

        for (std::vector<space_type>::iterator itr(spaces_.begin()); itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_particle(pid))
                return *itr;
        }
        throw NotFound("There is no particle having the given ParticleID.");
    }

    const_space_type get_space(const ParticleID& pid) const
    {
        if (root_->has_particle(pid))
            return const_space_type(root_, 0, 0);

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_particle(pid))
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

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            num_species += (*itr).num_species();
        }

        return num_species;
    }

    bool has_species(const Species& species) const
    {
        if (root_->has_species(species))
            return true;

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_species(species))
                return true;
        }
        return false;
    }

    // bool has_species_exact(const Species &sp) const;

    Integer num_molecules(const Species& species) const
    {
        Integer num_molecules(root_->num_molecules(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            num_molecules += (*itr).num_molecules(species);
        }
        return num_molecules;
    }

    Integer num_molecules_exact(const Species& species) const
    {
        Integer num_molecules(root_->num_molecules_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            num_molecules += (*itr).num_molecules_exact(species);
        }
        return num_molecules;
    }

    Real get_value(const Species& species) const
    {
        Real value(root_->get_value(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            value += (*itr).get_value(species);
        }
        return value;
    }

    Real get_value_exact(const Species& species) const
    {
        Real value(root_->get_value_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            value += (*itr).get_value_exact(species);
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
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            num_particles += (*itr).num_particles();
        }
        return num_particles;
    }

    Integer num_particles(const Species& species) const
    {
        Integer num_particles(root_->num_particles(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            num_particles += (*itr).num_particles(species);
        }
        return num_particles;
    }

    Integer num_particles_exact(const Species& species) const
    {
        Integer num_particles(root_->num_particles_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            num_particles += (*itr).num_particles_exact(species);
        }
        return num_particles;
    }

    bool has_particle(const ParticleID& pid) const
    {
        if (root_->has_particle(pid))
            return true;

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_particle(pid))
                return true;
        }
        return false;
    }

    identified_particle get_particle(const ParticleID& pid) const
    {
        if (root_->has_particle(pid))
            return root_->get_particle(pid);

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_particle(pid))
                return (*itr).get_particle(pid);
        }
        throw "Not Found";
    }

    std::vector<identified_particle> list_particles() const
    {
        std::vector<identified_particle> retval(root_->list_particles());
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<identified_particle> particles((*itr).list_particles());
            retval.insert(retval.end(), particles.begin(), particles.end());
        }
        return retval;
    }

    std::vector<identified_particle> list_particles(const Species& species) const
    {
        std::vector<identified_particle> retval(root_->list_particles(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<identified_particle> particles((*itr).list_particles(species));
            retval.insert(retval.end(), particles.begin(), particles.end());
        }
        return retval;
    }

    std::vector<identified_particle> list_particles_exact(const Species& species) const
    {
        std::vector<identified_particle> retval(root_->list_particles_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<identified_particle> particles((*itr).list_particles_exact(species));
            retval.insert(retval.end(), particles.begin(), particles.end());
        }
        return retval;
    }

    /*
     * LatticeSpaceTraits
     */
    const Real3 coordinate2position(const coordinate_type& coordinate) const
    {
        return get_space(coordinate).coordinate2position(coordinate);
    }

    coordinate_type position2coordinate(const Real3& position) const
    {
        return root_->position2coordinate(position);
    }

    Real get_volume(const Species& species) const
    {
        if (root_->has_molecule_pool(species))
        {
            if (root_->find_molecule_pool(species)->is_structure())
                return 0.0;
            else
                return root_->get_volume(species);
        }

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_molecule_pool(species))
            {
                if ((*itr).find_molecule_pool(species)->is_structure())
                    return 0.0;
                else
                    return (*itr).get_volume(species);
            }
        }

        return 0.0;
    }

    Integer num_neighbors(const coordinate_type& coordinate) const
    {
        return get_space(coordinate).num_neighbors(coordinate);
    }

    coordinate_type get_neighbor(coordinate_type coordinate, Integer nrand) const
    {
        return get_space(coordinate).get_neighbor(coordinate, nrand);
    }

    coordinate_type get_neighbor_boundary(coordinate_type coordinate, Integer nrand) const
    {
        return get_space(coordinate).get_neighbor_boundary(coordinate, nrand);
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
        return get_space_from_inner(coordinate).inner2coordinate(coordinate);
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
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<Species> new_species((*itr).list_species());
            species.insert(species.end(), new_species.begin(), new_species.end());
        }
        return species;
    }

    Integer num_voxels() const
    {
        Integer num_voxels(root_->num_voxels());
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            num_voxels += (*itr).num_voxels();
        }
        return num_voxels;
    }

    Integer num_voxels(const Species& species) const
    {
        Integer num_voxels(root_->num_voxels(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            num_voxels += (*itr).num_voxels(species);
        }
        return num_voxels;
    }

    Integer num_voxels_exact(const Species& species) const
    {
        Integer num_voxels(root_->num_voxels_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            num_voxels += (*itr).num_voxels_exact(species);
        }
        return num_voxels;
    }

    bool has_voxel(const ParticleID& pid) const
    {
        if (root_->has_voxel(pid))
            return true;

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_voxel(pid))
                return true;
        }
        return false;
    }

    identified_voxel get_voxel(const ParticleID& pid) const
    {
        return get_space(pid).get_voxel(pid);
    }

    std::vector<identified_voxel> list_voxels() const
    {
        std::vector<identified_voxel> voxels(root_->list_voxels());
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<identified_voxel> new_voxels((*itr).list_voxels());
            voxels.insert(voxels.end(), new_voxels.begin(), new_voxels.end());
        }
        return voxels;
    }

    std::vector<identified_voxel> list_voxels(const Species& species) const
    {
        std::vector<identified_voxel> voxels(root_->list_voxels(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<identified_voxel> new_voxels((*itr).list_voxels(species));
            voxels.insert(voxels.end(), new_voxels.begin(), new_voxels.end());
        }
        return voxels;
    }

    std::vector<identified_voxel> list_voxels_exact(const Species& species) const
    {
        std::vector<identified_voxel> voxels(root_->list_voxels_exact(species));
        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            std::vector<identified_voxel> new_voxels((*itr).list_voxels_exact(species));
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
            if ((*itr).has_species(species))
                return (*itr).find_voxel_pool(species);
        }
        throw NotFound("Not Found");
    }

    const VoxelPool* find_voxel_pool(const Species& species) const
    {
        if (root_->has_species(species))
            return root_->find_voxel_pool(species);

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_species(species))
                return (*itr).find_voxel_pool(species);
        }
        throw NotFound("Not Found");
    }

    bool has_molecule_pool(const Species& species) const
    {
        if (root_->has_molecule_pool(species))
            return true;

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_molecule_pool(species))
                return true;
        }
        return false;
    }

    std::pair<const MoleculePool*, coordinate_type>
    find_molecule_pool(const Species& species) const
    {
        if (root_->has_molecule_pool(species))
            return std::make_pair(root_->find_molecule_pool(species), 0);

        for (std::vector<space_type>::const_iterator itr(spaces_.begin());
             itr != spaces_.end(); ++itr)
        {
            if ((*itr).has_molecule_pool(species))
                return std::make_pair((*itr).find_molecule_pool(species), (*itr).offset());
        }
        throw "Not Found";
    }

    bool remove_particle(const ParticleID& pid)
    {
        space_type item(get_space_mut(pid));
        return item.remove_particle(pid);
    }

    bool on_structure(const Voxel& voxel)
    {
        return get_space_mut(voxel.coordinate()).on_structure(voxel);
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
        return get_space(coordinate).get_voxel_at(coordinate);
    }

    const VoxelPool* get_voxel_pool_at(const coordinate_type& coordinate) const
    {
        return get_space(coordinate).get_voxel_pool_at(coordinate);
    }

    /* Not wrapped function
    virtual const Particle particle_at(const coordinate_type& coord) const;
    */

    bool update_voxel(const ParticleID& pid, const Voxel& voxel)
    {
        // return get_space_mut(pid).update_voxel(pid, voxel);
        return get_space_mut(voxel.coordinate()).update_voxel(pid, voxel);
    }

    bool remove_voxel(const ParticleID pid)
    {
        return get_space_mut(pid).remove_voxel(pid);
    }

    bool remove_voxel(const coordinate_type& coordinate)
    {
        return get_space_mut(coordinate).remove_voxel(coordinate);
    }

    bool can_move(const coordinate_type& src, coordinate_type& dest)
    {
        boost::optional<const std::vector<coordinate_type>&> adjoinings(interfaces_.find(dest));

        if (adjoinings)
        {
            dest = pick(adjoinings.get(), rng());
            return false;
        }

        return get_space(src).can_move(src, dest);
    }

    bool move(const coordinate_type& src, const coordinate_type& dest,
              const std::size_t candidate=0)
    {
        return get_space_mut(src).move(src, dest, candidate);
    }

    /**
     * static members
     */
    static inline Real calculate_voxel_volume(const Real r)
    {
        return LatticeSpace::calculate_voxel_volume(r);
    }

    static inline Real3 calculate_hcp_lengths(const Real voxel_radius)
    {
        return LatticeSpace::calculate_hcp_lengths(voxel_radius);
    }

    static inline Integer3 calculate_shape(const Real3& edge_lengths, const Real& voxel_radius)
    {
        return LatticeSpace::calculate_shape(edge_lengths, voxel_radius, true);
    }

    static inline Real calculate_volume(const Real3& edge_lengths, const Real& voxel_radius)
    {
        return LatticeSpace::calculate_volume(edge_lengths, voxel_radius, true);
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

inline const Real
calculate_dimensional_factor(const VoxelPool* vp_a,
                             const VoxelPool* vp_b,
                             const boost::shared_ptr<const SpatiocyteWorld>& world)
{
    const Real voxel_radius(world->voxel_radius());
    const Real D_a(vp_a->D()),
               D_b(vp_b->D());
    const Shape::dimension_kind dim_a(vp_a->get_dimension()),
                                dim_b(vp_b->get_dimension());

    if (dim_a == Shape::THREE && dim_b == Shape::THREE)
    {
        return 2. / (3. * world->unit_voxel_volume() * (D_a + D_b) * voxel_radius);
    }
    else if (dim_a == Shape::TWO && dim_b == Shape::TWO)
    {
        const Real gamma(pow(2 * sqrt(2.0) + 4 * sqrt(3.0) + 3 * sqrt(6.0) + sqrt(22.0), 2) /
                         (72 * (6 * sqrt(2.0) + 4 * sqrt(3.0) + 3 * sqrt(6.0))));
        return gamma / (D_a + D_b);
    }
    else if (dim_a == Shape::THREE && dim_b == Shape::TWO)
    {
        const Real factor(sqrt(2.0) / (3 * D_a * voxel_radius));
        if (vp_b->is_structure())
            return factor * world->unit_area();
        return factor;
    }
    else if (dim_a == Shape::TWO && dim_b == Shape::THREE)
    {
        const Real factor(sqrt(2.0) / (3 * D_b * voxel_radius));
        if (vp_a->is_structure())
            return factor * world->unit_area();
        return factor;
    }
    throw NotSupported("The dimension of a structure must be two or three.");
}

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

#endif /* ECELL4_LATTICE_LATTICE_WORLD_HPP */
