#include "Context.hpp"
#include "MolecularType.hpp"
#include "StructureType.hpp"
#include "InterfaceType.hpp"
#include "LatticeSpaceVectorImpl.hpp"

namespace ecell4 {

LatticeSpaceVectorImpl::LatticeSpaceVectorImpl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const bool is_periodic) :
    base_type(edge_lengths, voxel_radius, is_periodic), is_periodic_(is_periodic)
{
    std::stringstream ss;
    ss << voxel_radius_;
    border_ = new MolecularType(Species("Border", ss.str(), "0"), vacant_);
    periodic_ = new MolecularType(Species("Periodic", ss.str(), "0"), vacant_);

    initialize_voxels(is_periodic_);
}

LatticeSpaceVectorImpl::~LatticeSpaceVectorImpl()
{
    delete border_;
    delete periodic_;
}

void LatticeSpaceVectorImpl::initialize_voxels(const bool is_periodic)
{
    const coordinate_type voxel_size(
        col_size_ * row_size_ * layer_size_);
    // std::cout << "voxel_size = " << voxel_size << std::endl;

    voxel_pools_.clear();
    molecule_pools_.clear();
    voxels_.clear();
    voxels_.reserve(voxel_size);
    for (coordinate_type coord(0); coord < voxel_size; ++coord)
    {
        if (!is_inside(coord))
        {
            if (is_periodic)
            {
                voxels_.push_back(periodic_);
            }
            else
            {
                voxels_.push_back(border_);
            }
        }
        else
        {
            voxels_.push_back(vacant_);
        }
    }
}

bool LatticeSpaceVectorImpl::update_structure(const Particle& p)
{
    //XXX: Particle does not have a location.
    Voxel v(p.species(), position2coordinate(p.position()), p.radius(), p.D());
    return update_voxel(ParticleID(), v);
}

/*
 * original methods
 */

const Species& LatticeSpaceVectorImpl::find_species(std::string name) const
{
    for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
         itr != voxel_pools_.end(); ++itr)
    {
        if ((*itr).first.serial() == name)
        {
            return (*itr).first;
        }
    }

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        if ((*itr).first.serial() == name)
        {
            return (*itr).first;
        }
    }
    throw NotFound(name);
}

std::vector<LatticeSpaceVectorImpl::coordinate_type>
LatticeSpaceVectorImpl::list_coords_exact(const Species& sp) const
{
    std::vector<coordinate_type> retval;

    molecule_pool_map_type::const_iterator itr(molecule_pools_.find(sp));
    if (itr == molecule_pools_.end())
    {
        return retval;
    }

    const boost::shared_ptr<MoleculePool>& vp((*itr).second);

    for (MoleculePool::const_iterator itr(vp->begin()); itr != vp->end(); ++itr)
    {
        retval.push_back((*itr).coordinate);
    }
    return retval;
}

std::vector<LatticeSpaceVectorImpl::coordinate_type>
LatticeSpaceVectorImpl::list_coords(const Species& sp) const
{
    std::vector<coordinate_type> retval;
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        if (!spmatch(sp, (*itr).first))
        {
            continue;
        }

        const boost::shared_ptr<MoleculePool>& vp((*itr).second);

        for (MoleculePool::const_iterator vitr(vp->begin());
             vitr != vp->end(); ++vitr)
        {
            retval.push_back((*vitr).coordinate);
        }
    }
    return retval;
}

/*
 * Protected functions
 */

LatticeSpaceVectorImpl::coordinate_type LatticeSpaceVectorImpl::get_coord(
    const ParticleID& pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        for (MoleculePool::const_iterator vitr(vp->begin());
             vitr != vp->end(); ++vitr)
        {
            if ((*vitr).pid == pid)
            {
                return (*vitr).coordinate;
            }
        }
    }
    return -1; //XXX: a bit dirty way
}

const VoxelPool* LatticeSpaceVectorImpl::get_voxel_pool_at(const coordinate_type& coord) const
{
    return voxels_.at(coord);
}

VoxelPool* LatticeSpaceVectorImpl::get_voxel_pool_at(const coordinate_type& coord)
{
    return voxels_.at(coord);
}

bool LatticeSpaceVectorImpl::remove_voxel(const ParticleID& pid)
{
    for (molecule_pool_map_type::iterator i(molecule_pools_.begin());
         i != molecule_pools_.end(); ++i)
    {
        const boost::shared_ptr<MoleculePool>& vp((*i).second);
        MoleculePool::const_iterator j(vp->find(pid));
        if (j != vp->end())
        {
            const coordinate_type coord((*j).coordinate);
            if (!vp->remove_voxel_if_exists(coord))
            {
                return false;
            }

            voxel_container::iterator itr(voxels_.begin() + coord);
            (*itr) = vp->location();
            vp->location()->add_voxel(
                coordinate_id_pair_type(ParticleID(), coord));
            return true;
        }
    }
    return false;
}

bool LatticeSpaceVectorImpl::remove_voxel(const coordinate_type& coord)
{
    voxel_container::iterator itr(voxels_.begin() + coord);
    VoxelPool* vp(*itr);
    if (vp->is_vacant())
    {
        return false;
    }
    if (vp->remove_voxel_if_exists(coord))
    {
        (*itr) = vp->location();
        vp->location()->add_voxel(
            coordinate_id_pair_type(ParticleID(), coord));
        return true;
    }
    return false;
}

bool LatticeSpaceVectorImpl::move(
    const coordinate_type& src, const coordinate_type& dest, const std::size_t candidate)
{
    return move_(src, dest, candidate).second;
}

bool LatticeSpaceVectorImpl::can_move(
    const coordinate_type& src, const coordinate_type& dest) const
{
    if (src == dest)
        return false;

    const VoxelPool* src_vp(voxels_.at(src));
    if (src_vp->is_vacant())
        return false;

    VoxelPool* dest_vp(voxels_.at(dest));

    if (dest_vp == border_)
        return false;

    if (dest_vp == periodic_)
        dest_vp = voxels_.at(apply_boundary_(dest));

    return (dest_vp == src_vp->location());
}

std::pair<LatticeSpaceVectorImpl::coordinate_type, bool>
    LatticeSpaceVectorImpl::move_(
        coordinate_type from, coordinate_type to,
        const std::size_t candidate)
{
    if (from == to)
    {
        // std::cerr << " from == to ";
        return std::pair<coordinate_type, bool>(from, false);
    }

    VoxelPool* from_vp(voxels_.at(from));
    if (from_vp->is_vacant())
    {
        return std::pair<coordinate_type, bool>(from, true);
    }

    VoxelPool* to_vp(voxels_.at(to));

    if (to_vp == border_)
    {
        // std::cerr << " to_vp is border ";
        return std::pair<coordinate_type, bool>(from, false);
    }
    else if (to_vp == periodic_)
    {
        to = apply_boundary_(to);
        to_vp = voxels_.at(to);
    }

    if (to_vp != from_vp->location())
    {
        // std::cerr << " to_vp is " << to_vp->species().serial() << " ";
        return std::pair<coordinate_type, bool>(to, false);
    }

    from_vp->replace_voxel(from, to, candidate);
    voxel_container::iterator from_itr(voxels_.begin() + from);
    (*from_itr) = to_vp;

    // to_vp->replace_voxel(to, coordinate_id_pair_type(ParticleID(), from));
    to_vp->replace_voxel(to, from);
    voxel_container::iterator to_itr(voxels_.begin() + to);
    (*to_itr) = from_vp;

    return std::pair<coordinate_type, bool>(to, true);
}

std::pair<LatticeSpaceVectorImpl::coordinate_type, bool>
    LatticeSpaceVectorImpl::move_(
        coordinate_id_pair_type& info, coordinate_type to)
{
    const coordinate_type from(info.coordinate);
    if (from == to)
    {
        return std::pair<coordinate_type, bool>(from, false);
    }

    VoxelPool* from_vp(voxels_.at(from));
    if (from_vp->is_vacant())
    {
        return std::pair<coordinate_type, bool>(from, true);
    }

    VoxelPool* to_vp(voxels_.at(to));

    if (to_vp == border_)
    {
        return std::pair<coordinate_type, bool>(from, false);
    }
    else if (to_vp == periodic_)
    {
        to = apply_boundary_(to);
        to_vp = voxels_.at(to);
    }

    if (to_vp != from_vp->location())
    {
        return std::pair<coordinate_type, bool>(to, false);
    }

    info.coordinate = to;
    voxel_container::iterator from_itr(voxels_.begin() + from);
    (*from_itr) = to_vp;

    // to_vp->replace_voxel(to, coordinate_id_pair_type(ParticleID(), from));
    to_vp->replace_voxel(to, from);
    voxel_container::iterator to_itr(voxels_.begin() + to);
    (*to_itr) = from_vp;

    return std::pair<coordinate_type, bool>(to, true);
}

const Particle LatticeSpaceVectorImpl::particle_at(
    const coordinate_type& coord) const
{
    const VoxelPool* vp(voxels_.at(coord));
    return Particle(
        vp->species(),
        coordinate2position(coord),
        vp->radius(), vp->D());
}

LatticeSpaceVectorImpl::identified_voxel
LatticeSpaceVectorImpl::get_voxel_at(const coordinate_type& coord) const
{
    const VoxelPool* vp(voxels_[coord]);
    const std::string loc(get_location_serial(vp));
    return std::make_pair(
        vp->get_particle_id(coord),
        Voxel(vp->species(), coord, vp->radius(), vp->D(), loc));
}

// /**
//  * Change the Species at v.coordinate() to v.species.
//  * The ParticleID must be kept after this update.
//  */
// void LatticeSpaceVectorImpl::update_voxel(const Voxel& v)
// {
//     const coordinate_type coord(v.coordinate());
//     VoxelPool* src_vp(voxels_.at(coord));
//     VoxelPool* new_vp(get_voxel_pool(v)); //XXX: need MoleculeInfo
// 
//     if (src_vp->with_voxels() != new_vp->with_voxels())
//     {
//         throw NotSupported("ParticleID is needed/lost.");
//     }
// 
//     new_vp->add_voxel(src_vp->pop(coord));
//     voxel_container::iterator itr(voxels_.begin() + coord);
//     (*itr) = new_vp;
// }

/*
 * Change the Species and coordinate of a Voxel with ParticleID, pid, to
 * v.species() and v.coordinate() respectively and return false.
 * If no Voxel with pid is found, create a new Voxel at v.coordiante() and return ture.
 */
bool LatticeSpaceVectorImpl::update_voxel(const ParticleID& pid, const Voxel& v)
{
    const LatticeSpaceVectorImpl::coordinate_type& to_coord(v.coordinate());
    if (!is_in_range(to_coord))
    {
        throw NotSupported("Out of bounds");
    }

    VoxelPool* new_vp(get_voxel_pool(v)); //XXX: need MoleculeInfo
    VoxelPool* dest_vp(get_voxel_pool_at(to_coord));

    if (dest_vp != new_vp->location())
    {
        throw NotSupported(
            "Mismatch in the location. Failed to place '"
            + new_vp->species().serial() + "' to '"
            + dest_vp->species().serial() + "'.");
    }

    const LatticeSpaceVectorImpl::coordinate_type
        from_coord(pid != ParticleID() ? get_coord(pid) : -1);
    if (from_coord != -1)
    {
        // move
        VoxelPool* src_vp(voxels_.at(from_coord));
        src_vp->remove_voxel_if_exists(from_coord);

        //XXX: use location?
        dest_vp->replace_voxel(to_coord, from_coord);
        voxel_container::iterator from_itr(voxels_.begin() + from_coord);
        (*from_itr) = dest_vp;

        new_vp->add_voxel(coordinate_id_pair_type(pid, to_coord));
        voxel_container::iterator to_itr(voxels_.begin() + to_coord);
        (*to_itr) = new_vp;
        return false;
    }

    // new
    dest_vp->remove_voxel_if_exists(to_coord);

    new_vp->add_voxel(coordinate_id_pair_type(pid, to_coord));
    voxel_container::iterator to_itr(voxels_.begin() + to_coord);
    (*to_itr) = new_vp;
    return true;
}

bool LatticeSpaceVectorImpl::add_voxels(const Species sp, std::vector<std::pair<ParticleID, coordinate_type> > voxels)
{
    // this function doesn't check location.
    VoxelPool *vp;
    try
    {
        vp = find_voxel_pool(sp);
    }
    catch (NotFound &e)
    {
        return false;
    }
    for (std::vector<std::pair<ParticleID, coordinate_type> >::iterator itr(voxels.begin());
            itr != voxels.end(); ++itr)
    {
        const ParticleID pid((*itr).first);
        const coordinate_type coord((*itr).second);
        VoxelPool* src_vp(get_voxel_pool_at(coord));
        src_vp->remove_voxel_if_exists(coord);
        vp->add_voxel(coordinate_id_pair_type(pid, coord));
        voxel_container::iterator vitr(voxels_.begin() + coord);
        (*vitr) = vp;
    }
    return true;
}

Integer LatticeSpaceVectorImpl::count_voxels(const boost::shared_ptr<VoxelPool>& vp) const
{
    return static_cast<Integer>(
        std::count(voxels_.begin(), voxels_.end(), vp.get()));
}

} // ecell4
