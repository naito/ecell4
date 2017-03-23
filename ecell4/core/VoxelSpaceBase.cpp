#include "VoxelSpaceBase.hpp"
#include "Context.hpp"
#include "InterfaceType.hpp"

namespace ecell4
{

VoxelSpaceBase::VoxelSpaceBase(const Real& voxel_radius)
    : t_(0.0),
      voxel_radius_(voxel_radius),
      vacant_(StructureType::allocVacant("Vacant", Shape::THREE))
{}

VoxelSpaceBase::~VoxelSpaceBase()
{
    delete vacant_;
}

/*
 * for LatticeSpace
 */
bool VoxelSpaceBase::on_structure(const Voxel& v)
{
    return get_voxel_pool_at(v.coordinate()) != get_voxel_pool(v)->location();
}

bool VoxelSpaceBase::make_structure_type(const Species& sp,
    Shape::dimension_kind dimension, const std::string loc)
{
    return make_location_type<StructureType>(sp, dimension, loc);
}

bool VoxelSpaceBase::make_interface_type(const Species& sp,
    Shape::dimension_kind dimension, const std::string loc)
{
    return make_location_type<InterfaceType>(sp, dimension, loc);
}

bool VoxelSpaceBase::make_molecular_pool(const Species& sp,
                                         Real radius,
                                         Real D,
                                         const std::string loc)
{
    if (has_molecule_pool(sp))
    {
        return false;
    }

    if (has_voxel_pool(sp))
    {
        throw IllegalState(
            "The given species is already assigned to the VoxelPool with no voxels.");
    }

    VoxelPool* location(get_vp_from_serial(loc));

    boost::shared_ptr<MoleculePool> vp(
            new MolecularType(sp, location, radius, D));
    std::pair<molecule_pool_map_type::iterator, bool> ins_result(
            molecule_pools_.insert(molecule_pool_map_type::value_type(sp, vp)));
    if (!ins_result.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return ins_result.second;
}


/*
 * VoxelSpaceBaseTraits
 */

std::vector<Species> VoxelSpaceBase::list_species() const
{
    std::vector<Species> keys;
    utils::retrieve_keys(voxel_pools_, keys);
    utils::retrieve_keys(molecule_pools_, keys);
    return keys;
}

bool VoxelSpaceBase::has_voxel(const ParticleID& pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        if (vp->find(pid) != vp->end())
            return true;
    }
    return false;
}

Integer VoxelSpaceBase::num_voxels() const
{
    Integer count(0);

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        count += count_voxels((*itr).second);
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        count += (*itr).second->size();
    }

    return count;
}

Integer VoxelSpaceBase::num_voxels(const Species& species) const
{
    Integer count(0);
    SpeciesExpressionMatcher spexp(species);

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        if (spexp.match((*itr).first))
            count += count_voxels((*itr).second);
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        if (spexp.match((*itr).first))
            count += (*itr).second->size();
    }
    return count;
}

Integer VoxelSpaceBase::num_voxels_exact(const Species& species) const
{
    {
        voxel_pool_map_type::const_iterator itr(voxel_pools_.find(species));
        if (itr != voxel_pools_.end())
        {
            return count_voxels((*itr).second);
        }
    }
    {
        molecule_pool_map_type::const_iterator itr(molecule_pools_.find(species));
        if (itr != molecule_pools_.end())
        {
            const boost::shared_ptr<MoleculePool>& vp((*itr).second);
            return vp->size();  // upcast
        }
    }
    return 0;
}

std::vector<VoxelSpaceBase::identified_voxel>
VoxelSpaceBase::list_voxels() const
{
    std::vector<identified_voxel> retval;

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        push_voxels(retval, vp);
    }

    return retval;
}

std::vector<VoxelSpaceBase::identified_voxel>
VoxelSpaceBase::list_voxels(const Species& species) const
{
    std::vector<identified_voxel> retval;
    SpeciesExpressionMatcher spexp(species);

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
            itr != molecule_pools_.end(); ++itr)
        if (spexp.match((*itr).first))
            push_voxels(retval, (*itr).second);
    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
            itr != voxel_pools_.end(); ++itr)
        if (spexp.match((*itr).first))
            push_voxels(retval, (*itr).second);

    return retval;
}

std::vector<VoxelSpaceBase::identified_voxel>
VoxelSpaceBase::list_voxels_exact(const Species& species) const
{
    std::vector<identified_voxel> retval;

    {
        molecule_pool_map_type::const_iterator itr(molecule_pools_.find(species));
        if (itr != molecule_pools_.end())
        {
            push_voxels(retval, (*itr).second);
            return retval;
        }
    }
    {
        voxel_pool_map_type::const_iterator itr(voxel_pools_.find(species));
        if (itr != voxel_pools_.end())
        {
            push_voxels(retval, (*itr).second);
            return retval;
        }
    }
    return retval;
}

VoxelSpaceBase::identified_voxel
VoxelSpaceBase::get_voxel(const ParticleID& pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        MoleculePool::container_type::const_iterator target_itr(vp->find(pid));
        if (target_itr != vp->end())
            return std::make_pair(pid,
                                  Voxel((*itr).first,
                                        (*target_itr).coordinate,
                                        vp->radius(),
                                        vp->D(),
                                        get_location_serial(vp)));
    }

    throw NotFound("No Voxel corresponding to a given ParticleID is found.");
}

bool VoxelSpaceBase::has_voxel_pool(const Species& sp) const
{
    return voxel_pools_.find(sp) != voxel_pools_.end();
}

const VoxelPool* VoxelSpaceBase::find_voxel_pool(const Species& species) const
{
    voxel_pool_map_type::const_iterator itr(voxel_pools_.find(species));
    if (itr != voxel_pools_.end())
    {
        return (*itr).second.get();
    }
    return find_molecule_pool(species);  // upcast
}

VoxelPool* VoxelSpaceBase::find_voxel_pool(const Species& species)
{
    voxel_pool_map_type::iterator itr(voxel_pools_.find(species));
    if (itr != voxel_pools_.end())
    {
        return (*itr).second.get();
    }
    return find_molecule_pool(species);  // upcast
}

bool VoxelSpaceBase::has_molecule_pool(const Species& species) const
{
    return molecule_pools_.find(species) != molecule_pools_.end();
}

const MoleculePool* VoxelSpaceBase::find_molecule_pool(const Species& species) const
{
    molecule_pool_map_type::const_iterator itr(molecule_pools_.find(species));
    if (itr != molecule_pools_.end())
    {
        return (*itr).second.get();  // upcast
    }
    throw NotFound("MoleculePool not found.");
}

MoleculePool* VoxelSpaceBase::find_molecule_pool(const Species& species)
{
    molecule_pool_map_type::iterator itr(molecule_pools_.find(species));
    if (itr != molecule_pools_.end())
    {
        return (*itr).second.get();  // upcast
    }
    throw NotFound("MoleculePool not found.");
}

bool VoxelSpaceBase::remove_particle(const ParticleID& pid)
{
    return remove_voxel(pid);
}

/*
 * SpaceTraits
 */
const Real VoxelSpaceBase::t() const
{
    return t_;
}

void VoxelSpaceBase::set_t(const Real& t)
{
    if (t < 0.0)
    {
        throw std::invalid_argument("the time must be positive.");
    }
    t_ = t;
}

/*
 * CompartmentSpaceTraits
 */

const Real VoxelSpaceBase::volume() const
{
    return inner_size() * voxel_volume();
}

Integer VoxelSpaceBase::num_species() const
{
    return voxel_pools_.size() + molecule_pools_.size();
}

bool VoxelSpaceBase::has_species(const Species& species) const
{
    return molecule_pools_.find(species) != molecule_pools_.end()
        || voxel_pools_.find(species)    != voxel_pools_.end();
}

Integer VoxelSpaceBase::num_molecules(const Species& species) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(species);

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        const Integer cnt(sexp.count((*itr).first));
        if (cnt > 0)
        {
            const boost::shared_ptr<VoxelPool>& vp((*itr).second);
            count += count_voxels(vp) * cnt;
        }
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const Integer cnt(sexp.count((*itr).first));
        if (cnt > 0)
        {
            const boost::shared_ptr<MoleculePool>& vp((*itr).second);
            count += vp->size() * cnt;
        }
    }

    return count;
}

Integer VoxelSpaceBase::num_molecules_exact(const Species& species) const
{
    return num_voxels_exact(species);
}

Real VoxelSpaceBase::get_value(const Species& species) const
{
    return static_cast<Real>(num_molecules(species));
}

Real VoxelSpaceBase::get_value_exact(const Species& species) const
{
    return static_cast<Real>(num_molecules_exact(species));
}

/*
 * ParticleSpaceTraits
 */
Integer VoxelSpaceBase::num_particles() const
{
    return num_voxels();
}

Integer VoxelSpaceBase::num_particles(const Species& species) const
{
    return num_voxels(species);
}

Integer VoxelSpaceBase::num_particles_exact(const Species& species) const
{
    return num_voxels_exact(species);
}

bool VoxelSpaceBase::has_particle(const ParticleID& pid) const
{
    return has_voxel(pid);
}

VoxelSpaceBase::identified_particle
VoxelSpaceBase::get_particle(const ParticleID& pid) const
{
    const Voxel v(get_voxel(pid).second);
    return std::make_pair(pid, Particle(
        v.species(), coordinate2position(v.coordinate()), v.radius(), v.D()));
}

std::vector<VoxelSpaceBase::identified_particle>
VoxelSpaceBase::list_particles() const
{
    const std::vector<identified_voxel> voxels(list_voxels());

    std::vector<identified_particle> retval;
    retval.reserve(voxels.size());
    for (std::vector<identified_voxel>::const_iterator
        i(voxels.begin()); i != voxels.end(); ++i)
    {
        const ParticleID& pid((*i).first);
        const Particle p(particle_at((*i).second.coordinate()));
        retval.push_back(std::make_pair(pid, p));
    }
    return retval;
}

std::vector<VoxelSpaceBase::identified_particle>
VoxelSpaceBase::list_particles(const Species& species) const
{
    const std::vector<identified_voxel> voxels(list_voxels(species));

    std::vector<identified_particle> retval;
    retval.reserve(voxels.size());
    for (std::vector<identified_voxel>::const_iterator
        i(voxels.begin()); i != voxels.end(); ++i)
    {
        const ParticleID& pid((*i).first);
        const Particle p(particle_at((*i).second.coordinate()));
        retval.push_back(std::make_pair(pid, p));
    }
    return retval;
}

std::vector<VoxelSpaceBase::identified_particle>
VoxelSpaceBase::list_particles_exact(const Species& species) const
{
    const std::vector<identified_voxel>
        voxels(list_voxels_exact(species));

    std::vector<identified_particle> retval;
    retval.reserve(voxels.size());
    for (std::vector<identified_voxel>::const_iterator
        i(voxels.begin()); i != voxels.end(); ++i)
    {
        const ParticleID& pid((*i).first);
        const Particle p(particle_at((*i).second.coordinate()));
        retval.push_back(std::make_pair(pid, p));
    }
    return retval;
}

/*
 * protected functions
 */
VoxelPool* VoxelSpaceBase::get_voxel_pool(const Voxel& v)
{
    const Species& sp(v.species());

    {
        voxel_pool_map_type::iterator itr(voxel_pools_.find(sp));
        if (itr != voxel_pools_.end())
        {
            return (*itr).second.get();
        }
    }

    {
        molecule_pool_map_type::iterator itr(molecule_pools_.find(sp));
        if (itr != molecule_pools_.end())
        {
            return (*itr).second.get();  // upcast
        }
    }

    // Create a new molecular pool

    const bool suc = make_molecular_pool(sp, v.radius(), v.D(), v.loc());
    if (!suc)
    {
        throw IllegalState("never reach here");
    }

    molecule_pool_map_type::iterator i = molecule_pools_.find(sp);
    if (i == molecule_pools_.end())
    {
        throw IllegalState("never reach here");
    }
    return (*i).second.get();  // upcast
}


/*
 * private functions
 */
std::string
VoxelSpaceBase::get_location_serial(const boost::shared_ptr<VoxelPool>& voxel_pool) const
{
    return voxel_pool->location()->is_vacant() ?
        "" : voxel_pool->location()->species().serial();
}

void VoxelSpaceBase::push_voxels(std::vector<identified_voxel>& voxels,
                                 const boost::shared_ptr<MoleculePool>& voxel_pool) const
{
    const std::string location_serial(get_location_serial(voxel_pool));
    for (MoleculePool::const_iterator itr(voxel_pool->begin()); itr != voxel_pool->end(); ++itr)
        voxels.push_back(std::make_pair((*itr).pid,
                                         Voxel(voxel_pool->species(),
                                               (*itr).coordinate,
                                               voxel_pool->radius(),
                                               voxel_pool->D(),
                                               location_serial)));
}

void VoxelSpaceBase::push_voxels(std::vector<identified_voxel>& voxels,
                                 const boost::shared_ptr<VoxelPool>& voxel_pool) const
{
    const std::string location_serial(get_location_serial(voxel_pool));

    for (coordinate_type coord(0); coord < size(); coord++)
    {
        if (get_voxel_pool_at(coord) != voxel_pool.get())
            continue;

        voxels.push_back(get_voxel_at(coord));
    }
}

VoxelPool* VoxelSpaceBase::get_vp_from_serial(const std::string& serial)
{
    if (serial == "") return vacant_;

    const Species sp(serial);
    try
    {
        return find_voxel_pool(sp);
    }
    catch (const NotFound& err)
    {
        // XXX: A VoxelPool for the structure (location) must be allocated
        // XXX: before the allocation of a Species on the structure.
        // XXX: The VoxelPool cannot be automatically allocated at the time
        // XXX: because its MoleculeInfo is unknown.
        // XXX: LatticeSpaceVectorImpl::load will raise a problem about this issue.
        // XXX: In this implementation, the VoxelPool for a structure is
        // XXX: created with default arguments.
        boost::shared_ptr<MoleculePool> mt(
                new MolecularType(sp, vacant_, voxel_radius_, 0));
        std::pair<molecule_pool_map_type::iterator, bool> ins_result(
                molecule_pools_.insert(molecule_pool_map_type::value_type(sp, mt)));
        if (!ins_result.second)
        {
            throw AlreadyExists("never reach here. find_voxel_pool seems wrong.");
        }
        return (*ins_result.first).second.get();
    }
}

template<typename T>
bool VoxelSpaceBase::make_location_type(const Species& sp,
                                        Shape::dimension_kind dimension,
                                        const std::string loc)
{
    if (has_voxel_pool(sp))
    {
        return false;
    }

    if (has_molecule_pool(sp))
    {
        throw IllegalState(
                "The given species is already assigned to the MoleculePool.");
    }

    VoxelPool* location(get_vp_from_serial(loc));

    boost::shared_ptr<VoxelPool> vp(
            new T(sp, location, voxel_radius_, dimension));
    std::pair<voxel_pool_map_type::iterator, bool> ins_result(
            voxel_pools_.insert(voxel_pool_map_type::value_type(sp, vp)));

    if (!ins_result.second)
        throw AlreadyExists("never reach here.");

    return ins_result.second;
}

} // ecell4
