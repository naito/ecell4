#include <stdexcept>
#include <fstream>
#include <algorithm>

#include "SpatiocyteWorld.hpp"
#include <ecell4/core/extras.hpp>
#include <ecell4/core/LatticeSpaceCellListImpl.hpp>

namespace ecell4
{

namespace spatiocyte
{

SpatiocyteWorld* create_spatiocyte_world_cell_list_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const Integer3& matrix_sizes,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return new SpatiocyteWorld(
        new LatticeSpaceCellListImpl(edge_lengths, voxel_radius, matrix_sizes),
        rng);
}

SpatiocyteWorld* create_spatiocyte_world_vector_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return new SpatiocyteWorld(new LatticeSpaceVectorImpl(edge_lengths, voxel_radius), rng);
}

SpatiocyteWorld* create_spatiocyte_world_offlattice_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    const Integer3 lattice_size(Integer(edge_lengths[0]/(2.0*voxel_radius)),
                                Integer(edge_lengths[1]/(2.0*voxel_radius)),
                                Integer(edge_lengths[2]/(2.0*voxel_radius)));
    return new SpatiocyteWorld(create_hcp_offlattice_space(voxel_radius, lattice_size), rng);
}


void SpatiocyteWorld::add_space(VoxelSpaceBase *space)
{
    for (std::size_t i(0); i< space->size(); ++i)
    {
        const Real3 position(space->coordinate2position(i));
        std::vector<coordinate_type> interfaces;

        for (std::size_t j(0); j < root_->size(); ++j)
        {
            if (length(root_->coordinate2position(j) - position) < voxel_radius() * 2)
            {
                interfaces_.add(j, i + size_);
                interfaces.push_back(j);
            }
        }

        for (std::vector<coordinate_type>::const_iterator itr(interfaces.begin());
             itr != interfaces.end(); ++itr)
        {
            for (Integer k(0); k < root_->num_neighbors(*itr); ++k)
            {
                const coordinate_type neighbor(root_->get_neighbor(*itr, k));

                if (std::find(interfaces.begin(), interfaces.end(), neighbor) != interfaces.end())
                {
                    neighbors_.add(i, neighbor);
                }
            }
        }
    }

    spaces_.push_back(space_type(space, size_, inner_size_));

    size_ += space->size();
    inner_size_ += space->inner_size();
}

MoleculeInfo SpatiocyteWorld::get_molecule_info(const Species& sp) const
{
    const bool with_D(sp.has_attribute("D"));
    const bool with_radius(sp.has_attribute("radius"));
    const bool with_loc(sp.has_attribute("location"));

    Real radius(with_radius ? std::atof(sp.get_attribute("radius").c_str())
                            : voxel_radius());
    Real D(with_D ? std::atof(sp.get_attribute("D").c_str())
                  : 0.0);
    std::string loc(with_loc ? sp.get_attribute("location")
                             : "");

    if (!with_D || !with_radius)
    {
        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            Species attributed(bound_model->apply_species_attributes(sp));

            if (!with_D && attributed.has_attribute("D"))
            {
                D = std::atof(attributed.get_attribute("D").c_str());
            }

            if (!with_radius && attributed.has_attribute("radius"))
            {
                radius = std::atof(
                    attributed.get_attribute("radius").c_str());
            }

            if (!with_loc && attributed.has_attribute("location"))
            {
                loc = attributed.get_attribute("location");
            }
        }
    }

    MoleculeInfo info = {radius, D, loc};
    return info;
}

// XXX: Never called
Integer SpatiocyteWorld::add_neighbors(const Species& sp,
    const SpatiocyteWorld::coordinate_type center)
{
    Integer count(0);
    const MoleculeInfo info(get_molecule_info(sp));
    for (Integer i(0); i < 12; ++i)
    {
        const coordinate_type n(get_neighbor(center, i));
        if (new_voxel(Voxel(sp, n, info.radius, info.D, info.loc)).second)
            ++count;
        else
            throw "Error in add_neighbors()";
    }
    return count;

    // Integer count(0);
    // const MoleculeInfo info(get_molecule_info(sp));
    // std::vector<SpatiocyteWorld::coordinate_type> neighbors(
    //         get_neighbors(center));
    // for (std::vector<SpatiocyteWorld::coordinate_type>::iterator itr(
    //             neighbors.begin()); itr != neighbors.end(); itr++)
    // {
    //     if (new_voxel(Voxel(sp, *itr, info.radius, info.D, info.loc)).second)
    //         ++count;
    //     else
    //         throw "Error in add_neighbors()";
    // }
    // return count;
}

std::pair<SpatiocyteWorld::coordinate_type, bool>
SpatiocyteWorld::check_neighbor(const coordinate_type coord,
                                const std::string& loc)
{
    std::vector<coordinate_type> tmp;
    std::vector<coordinate_type> additional_neighbors(neighbors_.get(coord));

    tmp.reserve(12 + additional_neighbors.size());
    for (unsigned int rnd(0); rnd < 12; ++rnd)
    {
        const coordinate_type neighbor(get_neighbor(coord, rnd));
        const VoxelPool* mt(get_voxel_pool_at(neighbor));
        const std::string serial(mt->is_vacant() ? "" : mt->species().serial());
        if (serial == loc)
        {
            tmp.push_back(neighbor);
        }
    }

    for (std::vector<coordinate_type>::const_iterator itr(additional_neighbors.begin());
         itr != additional_neighbors.end(); ++itr)
    {
        const VoxelPool* mt(get_voxel_pool_at(*itr));
        const std::string serial(mt->is_vacant() ? "" : mt->species().serial());
        if (serial == loc)
        {
            tmp.push_back(*itr);
        }
    }

    if (tmp.size() == 0)
    {
        return std::make_pair(coord, false);
    }

    return std::make_pair(tmp[rng()->uniform_int(0, tmp.size() - 1)], true);

    // const Integer rnd(rng()->uniform_int(0, 11));
    // const coordinate_type neighbor(get_neighbor(coord, rnd));
    // bool flg = find_voxel_pool(neighbor)->is_vacant(); //XXX: loc
    // return std::make_pair(neighbor, flg);
}

/*
 * Python API
 */

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

std::pair<SpatiocyteWorld::identified_particle, bool>
SpatiocyteWorld::new_particle(const Particle& p)
{
    // ParticleID pid(sidgen_());
    // const bool is_succeeded(update_particle(pid, p));
    // return std::make_pair(get_particle(pid), is_succeeded);
    const MoleculeInfo minfo(get_molecule_info(p.species()));
    const Voxel v(p.species(), position2coordinate(p.position()), p.radius(), p.D(), minfo.loc);
    if (on_structure(v))
    {
        return std::make_pair(std::make_pair(ParticleID(), p), false);
    }
    const std::pair<identified_voxel, bool> retval = new_voxel(v);
    return std::make_pair(std::make_pair(retval.first.first, p), retval.second);
}

bool SpatiocyteWorld::update_particle(const ParticleID& pid, const Particle& p)
{
    const MoleculeInfo minfo(get_molecule_info(p.species()));
    return update_voxel(pid, Voxel(p.species(),
        position2coordinate(p.position()), p.radius(), p.D(), minfo.loc));
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
    const std::vector<Species> non_structure_species(
            list_non_structure_species());

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

std::pair<SpatiocyteWorld::identified_voxel, bool>
SpatiocyteWorld::new_voxel_structure(const Species& sp, const coordinate_type& coord)
{
    const MoleculeInfo minfo(get_molecule_info(sp));
    return new_voxel_structure(Voxel(sp, coord, minfo.radius, minfo.D, minfo.loc));
}

std::pair<SpatiocyteWorld::identified_voxel, bool>
SpatiocyteWorld::new_voxel_structure(const Voxel& v)
{
    const bool is_succeeded(update_voxel(ParticleID(), v));
    return std::make_pair(std::make_pair(ParticleID(), v), is_succeeded);
}

std::pair<SpatiocyteWorld::identified_voxel, bool>
SpatiocyteWorld::new_voxel_interface(const Species& sp, const coordinate_type& coord)
{
    const MoleculeInfo minfo(get_molecule_info(sp));
    const Voxel voxel(sp, coord, minfo.radius, minfo.D, minfo.loc);
    const ParticleID pid;
    return std::make_pair(std::make_pair(pid, voxel), update_voxel(pid, voxel));
}

Integer SpatiocyteWorld::add_structure(const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    const MoleculeInfo info(get_molecule_info(sp));
    root_->make_structure_type(sp, shape->dimension(), info.loc);

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

Integer
SpatiocyteWorld::add_structure3(const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    const MoleculeInfo info(get_molecule_info(sp));
    Integer count(0);
    for (coordinate_type inner(0); inner < inner_size(); ++inner) {
        const coordinate_type coord(inner2coordinate(inner));
        const Real L(shape->is_inside(coordinate2position(coord)));
        if (L > 0)
            continue;

        const Voxel v(sp, coord, info.radius, info.D, info.loc);
        if (new_voxel_structure(v).second)
            ++count;
    }
    return count;
}

Integer
SpatiocyteWorld::add_structure2(const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    const MoleculeInfo info(get_molecule_info(sp));
    Integer count(0);
    for (coordinate_type inner(0); inner < inner_size(); ++inner) {
        const coordinate_type coord(inner2coordinate(inner));
        if (!is_surface_voxel(coord, shape))
            continue;

        const Voxel v(sp, coord, info.radius, info.D, info.loc);
        if (new_voxel_structure(v).second)
            ++count;
    }
    return count;
}

bool SpatiocyteWorld::is_surface_voxel(const coordinate_type coord,
                                       const boost::shared_ptr<const Shape> shape) const
{
    const Real L(shape->is_inside(coordinate2position(coord)));
    if (L > 0 || L < -2 * voxel_radius())
        return false;

    for (Integer i(0); i < 12; ++i)
        if (shape->is_inside(coordinate2position(get_neighbor(coord, i))) > 0)
            return true;

    return false;
}


Integer SpatiocyteWorld::add_interface(const Species& sp)
{
    return root_->make_interface_type(sp, Shape::UNDEF, get_molecule_info(sp).loc);
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
        const coordinate_type coord(
                inner2coordinate(rng()->uniform_int(0, inner_size() - 1)));
        //XXX: just for consistency. rather use below
        // const coordinate_type coord(rng()->uniform_int(0, size() - 1));

        const Voxel v(sp, coord, info.radius, info.D, info.loc);

        if (on_structure(v))
        {
            continue;
        }
        else if (new_voxel(v).second)
        {
            ++count;
        }
    }
    return true;
}

bool SpatiocyteWorld::add_molecules(const Species& sp, const Integer& num,
                                    const boost::shared_ptr<const Shape> shape)
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
        const Voxel v(sp,
                      position2coordinate(pos),
                      info.radius,
                      info.D,
                      info.loc);

        if (on_structure(v))
        {
            continue;
        }
        else if (new_voxel(v).second)
        {
            ++count;
        }
    }
    return true;
}

void SpatiocyteWorld::remove_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    const MoleculePool* mtype(find_molecule_pool(sp));
    if (mtype->size() < num)
    {
        throw std::invalid_argument(
            "The number of molecules cannot be negative.");
    }

    Integer count(0);
    while (count < num)
    {
        const Integer idx(rng_->uniform_int(0, mtype->size() - 1));
        if (remove_voxel(mtype->at(idx).coordinate))
        {
            ++count;
        }
    }
}

void SpatiocyteWorld::bind_to(boost::shared_ptr<Model> model)
{
    if (boost::shared_ptr<Model> bound_model = lock_model())
    {
        if (bound_model.get() != model.get())
        {
            std::cerr << "Warning: Model already bound to SpatiocyteWorld"
                      << std::endl;
        }
    }

    model_ = model;
}


/*
 * Wrapped functions
 */

void SpatiocyteWorld::save(const std::string& filename) const
{
#ifdef WITH_HDF5
    boost::scoped_ptr<H5::H5File>
        fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
    rng_->save(fout.get());
    sidgen_.save(fout.get());
    boost::scoped_ptr<H5::Group>
        group(new H5::Group(fout->createGroup("LatticeSpace")));
    root_->save_hdf5(group.get());
    extras::save_version_information(fout.get(), "ecell4-spatiocyte-0.0-1");
#else
    throw NotSupported("This method requires HDF5. The HDF5 support is turned off.");
#endif
}

void SpatiocyteWorld::load(const std::string& filename)
{
#ifdef WITH_HDF5
    boost::scoped_ptr<H5::H5File>
        fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
    const H5::Group group(fin->openGroup("LatticeSpace"));
    root_->load_hdf5(group);
    sidgen_.load(*fin);
    rng_->load(*fin);
#else
    throw NotSupported("This method requires HDF5. The HDF5 support is turned off.");
#endif
}

} // spatiocyte

} // ecell4
