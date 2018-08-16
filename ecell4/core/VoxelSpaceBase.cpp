#include "VoxelSpaceBase.hpp"
#include "Context.hpp"

namespace ecell4
{

/*
 * CompartmentSpace Traits
 */

std::vector<Species> VoxelSpaceBase::list_species() const
{
    std::vector<Species> keys;
    utils::retrieve_keys(voxel_pools_, keys);
    utils::retrieve_keys(molecule_pools_, keys);
    return keys;
}

/*
 * VoxelSpace Traits
 */

bool VoxelSpaceBase::has_voxel(const ParticleID& pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const MoleculePool* vp((*itr).second);
        if (vp->find(pid) != vp->end())
            return true;
    }
    return false;
}

Integer VoxelSpaceBase::num_voxels_exact(const Species& sp) const
{
    {
        voxel_pool_map_type::const_iterator itr(voxel_pools_.find(sp));
        if (itr != voxel_pools_.end())
        {
            return itr->second->size();
        }
    }

    {
        molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
        if (itr != molecule_pools_.end())
        {
            const MoleculePool* vp((*itr).second);
            return itr->second->size();
        }
    }

    return 0;
}

Integer VoxelSpaceBase::num_voxels(const Species& sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        if (sexp.match((*itr).first))
        {
            count += itr->second->size();
        }
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        if (sexp.match((*itr).first))
        {
            const MoleculePool* vp((*itr).second);
            count += vp->size();
        }
    }
    return count;
}

Integer VoxelSpaceBase::num_voxels() const
{
    Integer count(0);

    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        count += itr->second->size();
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const MoleculePool* vp((*itr).second);
        count += vp->size();
    }

    return count;
}

void VoxelSpaceBase::push_voxels(std::vector<std::pair<ParticleID, ParticleVoxel> >& voxels,
        const MoleculePool* voxel_pool,
        const Species& species) const
{
    const std::string location_serial(get_location_serial(voxel_pool));
    for (MoleculePool::const_iterator i(voxel_pool->begin()); i != voxel_pool->end(); ++i)
        voxels.push_back(
                std::make_pair(
                    (*i).pid,
                    ParticleVoxel(species, (*i).coordinate, voxel_pool->radius(),
                        voxel_pool->D(), location_serial)));
}

std::vector<std::pair<ParticleID, ParticleVoxel> >
VoxelSpaceBase::list_voxels() const
{
    std::vector<std::pair<ParticleID, ParticleVoxel> > retval;

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const MoleculePool* vp((*itr).second);
        push_voxels(retval, vp, vp->species());
    }

    return retval;
}

std::vector<std::pair<ParticleID, ParticleVoxel> >
VoxelSpaceBase::list_voxels(const Species& sp) const
{
    std::vector<std::pair<ParticleID, ParticleVoxel> > retval;
    SpeciesExpressionMatcher sexp(sp);

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
            itr != molecule_pools_.end(); ++itr)
        if (sexp.match((*itr).first))
            push_voxels(retval, (*itr).second, sp);

    return retval;
}

std::vector<std::pair<ParticleID, ParticleVoxel> >
VoxelSpaceBase::list_voxels_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, ParticleVoxel> > retval;

    molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
        push_voxels(retval, (*itr).second, sp);
    return retval;
}

boost::optional<ParticleVoxel>
VoxelSpaceBase::find_voxel(const ParticleID& pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const Species& species((*itr).first);
        const MoleculePool* pool((*itr).second);

        MoleculePool::container_type::const_iterator j(pool->find(pid));
        if (j != pool->end())
        {
            return ParticleVoxel(species,
                                 (*j).coordinate,
                                 pool->radius(),
                                 pool->D(),
                                 get_location_serial(pool));
        }
    }

    return boost::none;
}

VoxelPool* VoxelSpaceBase::find_voxel_pool(const Species& sp)
{
    if (sp.serial() == "")
        return &vacant_;

    voxel_pool_map_type::iterator itr(voxel_pools_.find(sp));
    if (itr != voxel_pools_.end())
    {
        return (*itr).second;
    }
    return find_molecule_pool(sp);  // upcast
}

const VoxelPool* VoxelSpaceBase::find_voxel_pool(const Species& sp) const
{
    if (sp.serial() == "")
        return &vacant_;

    voxel_pool_map_type::const_iterator itr(voxel_pools_.find(sp));
    if (itr != voxel_pools_.end())
    {
        return (*itr).second;
    }
    return find_molecule_pool(sp);  // upcast
}

bool VoxelSpaceBase::has_molecule_pool(const Species& sp) const
{
    return (molecule_pools_.find(sp) != molecule_pools_.end());
}

MoleculePool* VoxelSpaceBase::find_molecule_pool(const Species& sp)
{
    molecule_pool_map_type::iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
    {
        return (*itr).second;  // upcast
    }
    throw NotFound("MoleculePool not found.");
}

const MoleculePool* VoxelSpaceBase::find_molecule_pool(const Species& sp) const
{
    molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
    {
        return (*itr).second;  // upcast
    }
    throw NotFound("MoleculePool not found.");
}

bool VoxelSpaceBase::make_molecular_type(
        const Species& sp, Real radius, Real D, const std::string loc)
{
    molecule_pool_map_type::iterator itr(molecule_pools_.find(sp));
    if (itr != molecule_pools_.end())
    {
        return false;
    }
    else if (voxel_pools_.find(sp) != voxel_pools_.end())
    {
        throw IllegalState(
            "The given species is already assigned to the VoxelPool with no voxels.");
    }

    MoleculePool* vp(new MoleculePool(sp, find_voxel_pool(Species(loc)), radius, D));

    std::pair<molecule_pool_map_type::iterator, bool>
        retval(molecule_pools_.insert(molecule_pool_map_type::value_type(sp, vp)));

    if (!retval.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return retval.second;
}

bool VoxelSpaceBase::make_structure_type(const Species& sp,
    Shape::dimension_kind dimension, const std::string loc)
{
    voxel_pool_map_type::iterator itr(voxel_pools_.find(sp));
    if (itr != voxel_pools_.end())
    {
        return false;
    }
    else if (molecule_pools_.find(sp) != molecule_pools_.end())
    {
        throw IllegalState(
            "The given species is already assigned to the MoleculePool.");
    }

    VoxelPool*
        vp(new StructureType(sp, find_voxel_pool(Species(loc)), voxel_radius_, dimension));
    std::pair<voxel_pool_map_type::iterator, bool>
        retval(voxel_pools_.insert(voxel_pool_map_type::value_type(sp, vp)));
    if (!retval.second)
    {
        throw AlreadyExists("never reach here.");
    }
    return retval.second;
}

VoxelPool*
VoxelSpaceBase::get_voxel_pool(ParticleVoxel voxel)
{
    const Species& sp(voxel.species);

    try
    {
        return find_voxel_pool(sp);
    }
    catch (const NotFound &_e)
    {
        // Create a new molecular pool

        if (!make_molecular_type(sp, voxel.radius, voxel.D, voxel.loc))
        {
            throw IllegalState("never reach here");
        }

        molecule_pool_map_type::iterator itr = molecule_pools_.find(sp);

        if (itr == molecule_pools_.end())
        {
            throw IllegalState("never reach here");
        }

        return itr->second;
    }
}

} // ecell4
