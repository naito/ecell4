#include <algorithm>
#include "Context.hpp"
#include "VacantType.hpp"
#include "MoleculePool.hpp"
#include "OffLatticeSpace.hpp"

namespace ecell4 {

OffLatticeSpace::OffLatticeSpace(const Real& voxel_radius)
    : base_type(voxel_radius),
      voxels_(),
      positions_(),
      adjoinings_()
{}

OffLatticeSpace::OffLatticeSpace(const Real& voxel_radius,
                                 const position_container& positions,
                                 const coordinate_pair_list_type& adjoining_pairs,
                                 const Shape::dimension_kind& dimension)
    : base_type(voxel_radius, dimension),
      voxels_(),
      positions_(),
      adjoinings_()
{
    reset(positions, adjoining_pairs);
}

OffLatticeSpace::~OffLatticeSpace() {}

void OffLatticeSpace::reset(const position_container& positions,
                            const coordinate_pair_list_type& adjoining_pairs)
{
    voxels_.clear();
    positions_.clear();
    adjoinings_.clear();

    const std::size_t size(positions.size());

    voxels_.resize(size, &vacant_);
    positions_.resize(size);
    adjoinings_.resize(size);

    for (coordinate_type coord(0); coord < size; ++coord)
    {
        vacant_.add_voxel(ParticleID(), coord);
    }

    std::copy(positions.begin(), positions.end(), positions_.begin());

    for (coordinate_pair_list_type::const_iterator itr(adjoining_pairs.begin());
         itr != adjoining_pairs.end(); ++itr)
    {
        const coordinate_type coord0((*itr).first);
        const coordinate_type coord1((*itr).second);

        if (is_in_range(coord0) && is_in_range(coord1))
        {
            adjoinings_.at(coord0).push_back(coord1);
            adjoinings_.at(coord1).push_back(coord0);
        }
        else
        {
            throw IllegalState("A given pair is invalid.");
        }
    }
}

boost::optional<OffLatticeSpace::coordinate_type>
OffLatticeSpace::get_coord(const ParticleID& pid) const
{
    if (pid == ParticleID())
        return boost::none;

    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const MoleculePool* vp((*itr).second);
        for (MoleculePool::const_iterator vitr(vp->begin());
             vitr != vp->end(); ++vitr)
        {
            if ((*vitr).pid == pid)
            {
                return (*vitr).coordinate;
            }
        }
    }
    // throw NotFound("A corresponding particle is not found");
    return boost::none;
}


/*
 * public functions
 */

// Same as LatticeSpaceVectorImpl
bool OffLatticeSpace::update_voxel(const ParticleID& pid, ParticleVoxel v)
{
    const coordinate_type& to_coord(v.coordinate);
    if (!is_in_range(to_coord))
        throw NotSupported("Out of bounds");

    VoxelPool* new_vp(get_voxel_pool(v));
    VoxelPool* dest_vp(get_voxel_pool_at(to_coord));

    if (dest_vp != new_vp->location())
    {
        throw NotSupported(
            "Mismatch in the location. Failed to place '"
            + new_vp->species().serial() + "' to '"
            + dest_vp->species().serial() + "'.");
    }

    if (boost::optional<coordinate_type> from_coord = get_coord(pid))
    {
        // move
        voxels_.at(*from_coord)
               ->remove_voxel_if_exists(*from_coord);

        //XXX: use location?
        dest_vp->replace_voxel(to_coord, *from_coord);
        voxels_.at(*from_coord) = dest_vp;

        new_vp->add_voxel(pid, to_coord);
        voxels_.at(to_coord) = new_vp;

        return false;
    }

    // new
    dest_vp->remove_voxel_if_exists(to_coord);

    new_vp->add_voxel(pid, to_coord);
    voxels_.at(to_coord) = new_vp;

    return true;
}

bool OffLatticeSpace::add_voxel(
        const Species& species, const ParticleID& pid, const coordinate_type& coord)
{
    VoxelPool* vpool(find_voxel_pool(species));
    VoxelPool* location(get_voxel_pool_at(coord));

    if (vpool->location() != location)
        return false;

    location->remove_voxel_if_exists(coord);
    vpool->add_voxel(pid, coord);
    voxels_.at(coord) = vpool;

    return true;
}

// Same as LatticeSpaceVectorImpl
bool OffLatticeSpace::remove_voxel(const ParticleID& pid)
{
    for (molecule_pool_map_type::iterator i(molecule_pools_.begin());
         i != molecule_pools_.end(); ++i)
    {
        MoleculePool* vp((*i).second);
        MoleculePool::const_iterator j(vp->find(pid));
        if (j != vp->end())
        {
            const coordinate_type coord((*j).coordinate);
            if (!vp->remove_voxel_if_exists(coord))
            {
                return false;
            }

            voxels_.at(coord) = vp->location();

            vp->location()->add_voxel(ParticleID(), coord);
            return true;
        }
    }
    return false;
}

// Same as LatticeSpaceVectorImpl
bool OffLatticeSpace::remove_voxel(const coordinate_type& coord)
{
    VoxelPool* vp(voxels_.at(coord));
    if (vp->is_vacant())
    {
        return false;
    }
    if (vp->remove_voxel_if_exists(coord))
    {
        voxels_.at(coord) = vp->location();
        vp->location()->add_voxel(ParticleID(), coord);
        return true;
    }
    return false;
}

bool OffLatticeSpace::can_move(const coordinate_type& src, const coordinate_type& dest) const
{
    if (src == dest) return false;

    const VoxelPool* src_vp(voxels_.at(src));
    if (src_vp->is_vacant()) return false;

    VoxelPool* dest_vp(voxels_.at(dest));

    return (voxels_.at(dest) == src_vp->location());
}

bool OffLatticeSpace::move(const coordinate_type& src, const coordinate_type& dest,
        const std::size_t candidate)
{
    if (src == dest) return false;

    VoxelPool* src_vp(voxels_.at(src));
    if (src_vp->is_vacant()) return true;

    VoxelPool* dest_vp(voxels_.at(dest));
    if (dest_vp != src_vp->location()) return false;

    src_vp->replace_voxel(src, dest, candidate);
    voxels_.at(src) = dest_vp;

    dest_vp->replace_voxel(dest, src);
    voxels_.at(dest) = src_vp;

    return true;
}

OffLatticeSpace::coordinate_type
OffLatticeSpace::position2coordinate(const Real3& pos) const
{
    coordinate_type coordinate(0);
    Real shortest_length = length(positions_.at(0) - pos);

    for (coordinate_type coord(1); coord < size(); ++coord)
    {
        const Real len(length(positions_.at(coord) - pos));
        if (len < shortest_length)
        {
            coordinate = coord;
            shortest_length = len;
        }
    }

    return coordinate;
}

} // ecell4
