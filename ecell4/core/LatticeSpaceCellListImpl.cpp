#include "LatticeSpaceCellListImpl.hpp"
#include "StructureType.hpp"


namespace ecell4
{

Integer LatticeSpaceCellListImpl::num_molecules(const Species& sp) const
{
    Integer count(0);
    SpeciesExpressionMatcher sexp(sp);

    // for (voxel_pool_map_type::const_iterator itr(voxel_pools_.begin());
    //      itr != voxel_pools_.end(); ++itr)
    // {
    //     const Integer cnt(sexp.count((*itr).first));
    //     if (cnt > 0)
    //     {
    //         const boost::shared_ptr<VoxelPool>& vp((*itr).second);
    //         count += count_voxels(vp) * cnt;
    //     }
    // }

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

/*
 * Change the Species and coordinate of a Voxel with ParticleID, pid, to
 * v.species() and v.coordinate() respectively and return false.
 * If no Voxel with pid is found, create a new Voxel at v.coordiante() and return ture.
 */
bool LatticeSpaceCellListImpl::update_voxel(const ParticleID& pid, const Voxel& v)
{
    const coordinate_type& to_coord(v.coordinate());
    if (!is_in_range(to_coord))
    {
        throw NotSupported("Out of bounds");
    }

    VoxelPool* new_vp(get_voxel_pool(v)); //XXX: need MoleculeInfo
    VoxelPool* dest_vp(get_voxel_pool_at(to_coord));

    if (dest_vp != new_vp->location())
    {
        throw NotSupported("Mismatch in the location.");
    }

    if (pid != ParticleID())
    {
        const std::pair<VoxelPool*, coordinate_type>
            target(__get_coordinate(pid));
        const coordinate_type& from_coord(target.second);
        if (from_coord != -1)
        {
            // move
            VoxelPool* src_vp(target.first);
            src_vp->remove_voxel_if_exists(from_coord);

            //XXX: use location?
            dest_vp->replace_voxel(to_coord, from_coord);

            new_vp->add_voxel(coordinate_id_pair_type(pid, to_coord));

            if (!dest_vp->is_vacant())
            {
                update_matrix(from_coord, dest_vp);
                update_matrix(to_coord, new_vp);
            }
            else
            {
                update_matrix(from_coord, to_coord, new_vp);
            }
            return true;
        }
    }

    // new
    dest_vp->remove_voxel_if_exists(to_coord);

    new_vp->add_voxel(coordinate_id_pair_type(pid, to_coord));
    update_matrix(to_coord, new_vp);
    return true;
}

std::pair<VoxelPool*, LatticeSpaceCellListImpl::coordinate_type>
    LatticeSpaceCellListImpl::__get_coordinate(const ParticleID& pid)
{
    for (molecule_pool_map_type::iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        MoleculePool::container_type::const_iterator j(vp->find(pid));
        if (j != vp->end())
        {
            return std::pair<VoxelPool*, coordinate_type>(
                vp.get(), (*j).coordinate);
        }
    }
    return std::make_pair<VoxelPool*, coordinate_type>(NULL, -1); //XXX: a bit dirty way
}

std::pair<const VoxelPool*, LatticeSpaceCellListImpl::coordinate_type>
    LatticeSpaceCellListImpl::__get_coordinate(const ParticleID& pid) const
{
    for (molecule_pool_map_type::const_iterator itr(molecule_pools_.begin());
         itr != molecule_pools_.end(); ++itr)
    {
        const boost::shared_ptr<MoleculePool>& vp((*itr).second);
        MoleculePool::container_type::const_iterator j(vp->find(pid));
        if (j != vp->end())
        {
            return std::pair<const VoxelPool*, coordinate_type>(
                vp.get(), (*j).coordinate);
        }
    }
    return std::make_pair<const VoxelPool*, coordinate_type>(NULL, -1); //XXX: a bit dirty way
}

const VoxelPool* LatticeSpaceCellListImpl::get_voxel_pool_at(
    const LatticeSpaceCellListImpl::coordinate_type& coord) const
{
    /**
     XXX: This may not work
     */
    if (!is_in_range(coord))
    {
        throw NotSupported("Out of bounds");
    }

    if (!is_inside(coord))
    {
        if (is_periodic_)
        {
            return periodic_;
        }
        else
        {
            return border_;
        }
    }

    // for (spmap::iterator itr(spmap_.begin());
    //     itr != spmap_.end(); ++itr)
    // {
    //     VoxelPool& vp((*itr).second);
    //     if (vp.is_vacant())
    //     {
    //         continue;
    //     }

    //     VoxelPool::container_type::const_iterator j(vp.find(coord));
    //     if (j != vp.end())
    //     {
    //         return (&vp);
    //     }
    // }

    const cell_type& cell(matrix_[coordinate2index(coord)]);
    if (cell.size() == 0)
    {
        return vacant_;
    }

    cell_type::const_iterator i(find_from_cell(coord, cell));
    if (i != cell.end())
    {
        return (*i).first;
    }

    return vacant_;
}

VoxelPool* LatticeSpaceCellListImpl::get_voxel_pool_at(
    const LatticeSpaceCellListImpl::coordinate_type& coord)
{
    /**
     XXX: This may not work
     */
    if (!is_in_range(coord))
    {
        throw NotSupported("Out of bounds");
    }

    if (!is_inside(coord))
    {
        if (is_periodic_)
        {
            return periodic_;
        }
        else
        {
            return border_;
        }
    }

    // for (spmap::iterator itr(spmap_.begin());
    //     itr != spmap_.end(); ++itr)
    // {
    //     VoxelPool& vp((*itr).second);
    //     if (vp.is_vacant())
    //     {
    //         continue;
    //     }

    //     VoxelPool::container_type::const_iterator j(vp.find(coord));
    //     if (j != vp.end())
    //     {
    //         return (&vp);
    //     }
    // }

    const cell_type& cell(matrix_[coordinate2index(coord)]);
    if (cell.size() == 0)
    {
        return vacant_;
    }

    cell_type::const_iterator i(find_from_cell(coord, cell));
    if (i != cell.end())
    {
        return (*i).first;
    }

    return vacant_;
}

void LatticeSpaceCellListImpl::add_structure(const Species& sp,
    const boost::shared_ptr<const Shape>& s, const std::string loc)
{
    make_structure_type(sp, s->dimension(), loc);

    structure_container_type::const_iterator i(structures_.find(sp));
    if (i != structures_.end())
    {
        throw NotSupported("not supported yet.");
    }
    structures_.insert(std::make_pair(sp, s));
}

const boost::shared_ptr<const Shape>& LatticeSpaceCellListImpl::get_structure(
    const Species& sp) const
{
    structure_container_type::const_iterator i(structures_.find(sp));
    if (i == structures_.end())
    {
        throw NotFound("not found.");
    }
    return (*i).second;
}

const Shape::dimension_kind LatticeSpaceCellListImpl::get_structure_dimension(
    const Species& sp) const
{
    structure_container_type::const_iterator i(structures_.find(sp));
    if (i == structures_.end())
    {
        return Shape::THREE; // Default value (ex. for VACANT type)
    }
    return (*i).second->dimension();
}

} // ecell4
