#include <stdexcept>
#include <fstream>

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteWorld::add_space(space_type space)
{
    spaces_.push_back(space);

    for (std::size_t i(0); i < space->size(); ++i)
    {
        const Real3 position(space->coordinate2position(i));
        const Voxel nearest(get_root(), get_root()->position2coordinate(position));

        if (length(nearest.position() - position) < voxel_radius()*2)
        {
            interfaces_.add(nearest, Voxel(space, i));
        }
    }

    for (OneToManyMap<Voxel>::const_iterator itr(interfaces_.begin());
         itr != interfaces_.end(); ++itr)
    {
        const Voxel& interface((*itr).first);
        const std::vector<Voxel>& targets((*itr).second);

        std::vector<Voxel> neighbors;
        for (Integer i(0); i < interface.num_neighbors(); ++i)
        {
            const Voxel neighbor(interface.get_neighbor(i));
            if (! interfaces_.find(neighbor))
                neighbors.push_back(neighbor);
        }

        for (std::vector<Voxel>::const_iterator jtr(targets.begin());
             jtr != targets.end(); ++jtr)
        {
            neighbors_.extend(*jtr, neighbors);
        }
    }

    size_ += space->size();
}

void SpatiocyteWorld::set_value(const Species& sp, const Real value)
{
    const Integer num1 = static_cast<Integer>(value);
    const Integer num2 = num_molecules_exact(sp);
    if (num1 > num2)
    {
        add_molecules(sp, num1 - num2);
    }
    else if (num1 < num2)
    {
        remove_molecules(sp, num2 - num1);
    }
}

std::vector<std::pair<ParticleID, Particle> >
SpatiocyteWorld::list_structure_particles() const
{
    const std::vector<Species> structure_species(list_structure_species());

    typedef std::vector<std::vector<std::pair<ParticleID, Particle> > > tmp_type;
    tmp_type tmp_vector(structure_species.size());
    Integer num_elements;

    for (std::vector<Species>::const_iterator itr(structure_species.begin());
            itr != structure_species.end(); ++itr)
    {
        std::vector<std::pair<ParticleID, Particle> > tmp(list_particles(*itr));
        tmp_vector.push_back(tmp);
        num_elements += tmp.size();
    }

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(num_elements);
    for (tmp_type::const_iterator itr(tmp_vector.begin());
            itr != tmp_vector.end(); ++itr)
    {
        retval.insert(retval.end(), (*itr).begin(), (*itr).end());
    }

    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
SpatiocyteWorld::list_non_structure_particles() const
{
    const std::vector<Species> non_structure_species(list_non_structure_species());

    typedef std::vector<std::vector<std::pair<ParticleID, Particle> > > tmp_type;
    tmp_type tmp_vector(non_structure_species.size());
    Integer num_elements;

    for (std::vector<Species>::const_iterator itr(non_structure_species.begin());
            itr != non_structure_species.end(); ++itr)
    {
        std::vector<std::pair<ParticleID, Particle> > tmp(list_particles(*itr));
        tmp_vector.push_back(tmp);
        num_elements += tmp.size();
    }

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(num_elements);
    for (tmp_type::const_iterator itr(tmp_vector.begin());
            itr != tmp_vector.end(); ++itr)
    {
        retval.insert(retval.end(), (*itr).begin(), (*itr).end());
    }

    return retval;
}

std::vector<Species> SpatiocyteWorld::list_non_structure_species() const
{
    const std::vector<Species> species(list_species());
    std::vector<Species> retval;
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        if (!find_voxel_pool(*itr)->is_structure())
            retval.push_back(*itr);
    }
    return retval;
}

std::vector<Species> SpatiocyteWorld::list_structure_species() const
{
    const std::vector<Species> species(list_species());
    std::vector<Species> retval;
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        if (find_voxel_pool(*itr)->is_structure())
            retval.push_back(*itr);
    }
    return retval;
}

bool SpatiocyteWorld::add_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    const MoleculeInfo info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const Voxel voxel(coordinate2voxel(rng()->uniform_int(0, size()-1)));

        if (voxel.get_voxel_pool()->species().serial() != info.loc)
        {
            continue;
        }
        else if (new_voxel(sp, voxel))
        {
            ++count;
        }
    }
    return true;
}

bool SpatiocyteWorld::add_molecules(
    const Species& sp, const Integer& num, const boost::shared_ptr<const Shape> shape)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    const MoleculeInfo info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const Real3 pos(shape->draw_position(rng_));
        const Voxel voxel(position2voxel(pos));

        if (voxel.get_voxel_pool()->species().serial() != info.loc)
        {
            continue;
        }
        else if (new_voxel(sp, voxel))
        {
            ++count;
        }
    }
    return true;
}

Integer SpatiocyteWorld::add_structure(
    const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    const MoleculeInfo info(get_molecule_info(sp));
    get_root()->make_structure_type(sp, shape->dimension(), info.loc);

    switch (shape->dimension())
    {
    case Shape::THREE:
        return add_structure3(sp, shape);
    case Shape::TWO:
        return add_structure2(sp, shape);
    case Shape::ONE:
    case Shape::UNDEF:
        break;
    }

    throw NotSupported("The dimension of a shape must be two or three.");
}

Integer SpatiocyteWorld::add_structure3(const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    const MoleculeInfo info(get_molecule_info(sp));
    Integer count(0);
    for (coordinate_type coord(0); coord < size(); ++coord) {
        const Voxel voxel(coordinate2voxel(coord));
        const Real L(shape->is_inside(voxel.position()));
        if (L > 0)
            continue;

        if (voxel.get_voxel_pool()->species().serial() != info.loc)
        {
            continue;
        }

        if (new_voxel_structure(sp, voxel))
            ++count;
    }
    return count;
}

Integer
SpatiocyteWorld::add_structure2(
        const Species& sp,
        const boost::shared_ptr<const Shape> shape)
{
    const MoleculeInfo info(get_molecule_info(sp));
    Integer count(0);
    for (coordinate_type coord(0); coord < size(); ++coord) {
        const Voxel voxel(coordinate2voxel(coord));
        if (!is_surface_voxel(voxel, shape))
            continue;

        if (voxel.get_voxel_pool()->species().serial() != info.loc)
        {
            continue;
        }

        if (new_voxel_structure(sp, voxel))
            ++count;
    }
    return count;
}

bool
SpatiocyteWorld::is_surface_voxel(
        const Voxel& voxel,
        const boost::shared_ptr<const Shape> shape) const
{
    const Real L(shape->is_inside(voxel.position()));
    if (L > 0 || L < -2 * voxel_radius())
        return false;

    for (Integer i(0); i < voxel.num_neighbors(); ++i)
        if (shape->is_inside(voxel.get_neighbor(i).position()) > 0)
            return true;

    return false;
}

void SpatiocyteWorld::remove_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    boost::shared_ptr<const MoleculePool> mtype(find_molecule_pool(sp));
    if (mtype->size() < num)
    {
        throw std::invalid_argument(
            "The number of molecules cannot be negative.");
    }

    Integer count(0);
    while (count < num)
    {
        const Integer idx(rng_->uniform_int(0, mtype->size() - 1));
        if (coordinate2voxel(mtype->at(idx).coordinate).clear())
        {
            ++count;
        }
    }
}

boost::optional<Voxel>
SpatiocyteWorld::check_neighbor(const Voxel& voxel, const std::string& loc)
{
    const std::size_t num_neighbors(voxel.num_neighbors());

    std::vector<Voxel> tmp;
    tmp.reserve(num_neighbors);

    for (unsigned int rnd(0); rnd < num_neighbors; ++rnd)
    {
        const Voxel neighbor(voxel.get_neighbor(rnd));
        boost::shared_ptr<const VoxelPool> mt(neighbor.get_voxel_pool());
        const std::string serial(mt->is_vacant() ? "" : mt->species().serial());
        if (serial == loc)
        {
            tmp.push_back(neighbor);
        }
    }

    if (tmp.size() != 0)
    {
        return tmp[rng()->uniform_int(0, tmp.size()-1)];
    }

    if (boost::optional<const std::vector<Voxel>&> neighbors = find_neighbors(voxel))
    {
        tmp.reserve(neighbors->size());

        for (std::vector<Voxel>::const_iterator itr(neighbors->begin());
             itr != neighbors->end(); ++itr)
        {
            boost::shared_ptr<const VoxelPool> mt(itr->get_voxel_pool());
            const std::string serial(mt->is_vacant() ? "" : mt->species().serial());
            if (serial == loc)
            {
                tmp.push_back(*itr);
            }
        }

        if (tmp.size() != 0)
        {
            return tmp[rng()->uniform_int(0, tmp.size()-1)];
        }
    }

    return boost::none;
}

} // spatiocyte

} // ecell4
