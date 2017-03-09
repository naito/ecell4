#include "SpatiocyteEvent.hpp"
#include "utils.hpp"

namespace ecell4
{

namespace spatiocyte
{

inline ReactionInfo
apply_vanishment(boost::shared_ptr<SpatiocyteWorld> world,
                 const ReactionInfo::identified_voxel& p0,
                 const ReactionInfo::identified_voxel& p1)
{
    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(p0);
    rinfo.add_reactant(p1);

    world->remove_voxel(p0.second.coordinate());
    world->remove_voxel(p1.second.coordinate());

    return rinfo;
}

inline ReactionInfo
apply_ab2c(boost::shared_ptr<SpatiocyteWorld> world,
           const ReactionInfo::identified_voxel& p0,
           const ReactionInfo::identified_voxel& p1,
           const Species& product_species)
{
    // A and B (from_info and to_info) become C (product_species)
    const std::string& location(world->get_molecule_info(product_species).loc);
    const std::string& fserial(get_serial(world, p0.second.coordinate()));
    const std::string& floc(get_location(world, p0.second.coordinate()));
    const std::string& tserial(get_serial(world, p1.second.coordinate()));
    const std::string& tloc(get_location(world, p1.second.coordinate()));

    ReactionInfo rinfo(world->t());

    if (tserial == location || tloc == location)
    {
        // B is on the location of C, or the location itself.
        // Place C at the coordinate of B, and remove A.
        rinfo.add_reactant(p0);
        rinfo.add_reactant(p1);

        if (tserial != location)
        {
            world->remove_voxel(p1.second.coordinate());
        }

        world->remove_voxel(p0.second.coordinate());
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world->new_voxel(product_species, p1.second.coordinate()));

        rinfo.add_product(new_mol.first);
    }
    else if (fserial == location || floc == location)
    {
        // A is on the location of C, or the location itself.
        // Place C at the coordinate of A, and remove B.
        rinfo.add_reactant(p0);
        rinfo.add_reactant(p1);

        if (fserial != location)
        {
            world->remove_voxel(p0.second.coordinate());
        }

        world->remove_voxel(p1.second.coordinate());
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
            world->new_voxel(product_species, p0.second.coordinate()));

        rinfo.add_product(new_mol.first);
    }
    return rinfo;
}

inline ReactionInfo
apply_ab2cd_in_order(boost::shared_ptr<SpatiocyteWorld> world,
                     const ReactionInfo::identified_voxel& p0,
                     const ReactionInfo::identified_voxel& p1,
                     const Species& product_species0,
                     const Species& product_species1,
                     const SpatiocyteWorld::coordinate_type coord0,
                     const SpatiocyteWorld::coordinate_type coord1)
{
    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(p0);
    rinfo.add_reactant(p1);

    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
        world->new_voxel(product_species0, coord0));
    if (!new_mol0.second)
    {
        throw IllegalState("no place for " + product_species0.serial());
    }
    std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
        world->new_voxel(product_species1, coord1));
    if (!new_mol1.second)
    {
        throw IllegalState("no place for " + product_species1.serial());
    }

    rinfo.add_product(new_mol0.first);
    rinfo.add_product(new_mol1.first);

    return rinfo;
}

inline ReactionInfo
apply_ab2cd(boost::shared_ptr<SpatiocyteWorld> world,
            const ReactionInfo::identified_voxel& p0,
            const ReactionInfo::identified_voxel& p1,
            const Species& product_species0,
            const Species& product_species1)
{
    const SpatiocyteWorld::coordinate_type from_coord(p0.second.coordinate());
    const SpatiocyteWorld::coordinate_type to_coord(p1.second.coordinate());
    const std::string aserial(get_serial(world, from_coord));
    const std::string aloc(get_location(world, from_coord));
    const std::string bserial(get_serial(world, to_coord));
    const std::string bloc(get_location(world, to_coord));
    const std::string cloc(world->get_molecule_info(product_species0).loc);
    const std::string dloc(world->get_molecule_info(product_species1).loc);

    if (aserial == cloc || aloc == cloc)
    {
        if (bserial == dloc || bloc == dloc)
        {
            if (aserial != cloc)
            {
                // Remove A once if A is not the location of C
                world->remove_voxel(p0.second.coordinate());
            }
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                world, p0, p1, product_species0, product_species1,
                from_coord, to_coord);
        }
        else
        {
            std::pair<SpatiocyteWorld::coordinate_type, bool>
                neighbor(world->check_neighbor(to_coord, dloc));

            if (neighbor.second)
            {
                world->remove_voxel(p1.second.coordinate());
                if (aserial != cloc)
                {
                    // Remove A once if A is not the location of C
                    world->remove_voxel(p0.second.coordinate());
                }
                return apply_ab2cd_in_order(
                    world, p0, p1, product_species0, product_species1,
                    from_coord, neighbor.first);
            }
        }
    }
    else if (aserial == dloc || aloc == dloc)
    {
        if (bserial == cloc || bloc == dloc)
        {
            if (aserial != dloc)
            {
                // Remove A once if A is not the location of D
                world->remove_voxel(p0.second.coordinate());
            }
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                world, p0, p1, product_species0, product_species1,
                to_coord, from_coord);
        }
        else
        {
            std::pair<SpatiocyteWorld::coordinate_type, bool>
                neighbor(world->check_neighbor(to_coord, cloc));

            if (neighbor.second)
            {
                world->remove_voxel(p1.second.coordinate());
                if (aserial != dloc)
                {
                    // Remove A once if A is not the location of D
                    world->remove_voxel(p0.second.coordinate());
                }
                return apply_ab2cd_in_order(
                    world, p0, p1, product_species0, product_species1,
                    neighbor.first, from_coord);
            }
        }
    }
    else if (bserial == cloc || bloc == cloc)
    {
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(to_coord, dloc));

        if (neighbor.second)
        {
            world->remove_voxel(p0.second.coordinate());
            if (bserial != cloc)
            {
                // Remove B once if B is not the location of C
                world->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                world, p0, p1, product_species0, product_species1,
                to_coord, neighbor.first);
        }
    }
    else if (bserial == dloc || bloc == dloc)
    {
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(to_coord, dloc));

        if (neighbor.second)
        {
            world->remove_voxel(p0.second.coordinate());
            if (bserial != dloc)
            {
                // Remove B once if B is not the location of D
                world->remove_voxel(p1.second.coordinate());
            }
            return apply_ab2cd_in_order(
                world, p0, p1, product_species0, product_species1,
                neighbor.first, to_coord);
        }
    }

    return ReactionInfo(world->t());
}

inline ReactionInfo
apply_second_order_reaction(boost::shared_ptr<SpatiocyteWorld> world,
                            const ReactionRule& reaction_rule,
                            const ReactionInfo::identified_voxel& p0,
                            const ReactionInfo::identified_voxel& p1)
{
    const ReactionRule::product_container_type& products(reaction_rule.products());

    switch (products.size())
    {
        case 0:
            return apply_vanishment(world, p0, p1);
        case 1:
            return apply_ab2c(world, p0, p1, *(products.begin()));
        case 2:
            return apply_ab2cd(world, p0, p1,
                            *(products.begin()), *(++(products.begin())));
        default:
            return ReactionInfo(world->t());
    }
}

StepEvent::StepEvent(boost::shared_ptr<Model> model,
                     boost::shared_ptr<SpatiocyteWorld> world,
                     const Species& species, const Real& t, const Real alpha)
    : SpatiocyteEvent(t), model_(model), world_(world), species_(species), alpha_(alpha)
{
    const SpatiocyteWorld::molecule_info_type
        minfo(world_->get_molecule_info(species));
    const Real R(world_->voxel_radius());
    const Real D(minfo.D);
    const Real sqRperD(pow(R, 2.0)/D);
    const VoxelPool* mtype(world_->find_voxel_pool(species));
    if (D <= 0)
    {
        dt_ = inf;
    } else if(mtype->get_dimension() == Shape::THREE) {
        dt_ = 2 * sqRperD / 3 * alpha_;
    } else if(mtype->get_dimension() == Shape::TWO) {
        // TODO: Regular Lattice
        // dt_  = pow((2*sqrt(2.0)+4*sqrt(3.0)+3*sqrt(6.0)+sqrt(22.0))/
        //           (6*sqrt(2.0)+4*sqrt(3.0)+3*sqrt(6.0)), 2) * sqRperD * alpha_;
        dt_ = sqRperD * alpha_;
    } else if(mtype->get_dimension() == Shape::ONE) {
        dt_ = 2 * sqRperD * alpha_;
    }
    else
    {
        throw NotSupported(
            "The dimension of a structure must be two or three.");
    }

    time_ = t + dt_;
}

void StepEvent::fire_()
{
    walk(alpha_);
    time_ += dt_;
}

void StepEvent::walk(const Real& alpha)
{
    if (alpha < 0 || alpha > 1)
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    const MoleculePool* mtype(world_->find_molecule_pool(species_));

    if (mtype->get_dimension() == Shape::THREE)
        walk_in_space_(mtype, alpha);
    else // dimension == TWO, etc.
        walk_on_surface_(mtype, alpha);
}

void StepEvent::walk_in_space_(const MoleculePool* mtype, const Real& alpha)
{
    const boost::shared_ptr<RandomNumberGenerator> rng(world_->rng());
    MoleculePool::container_type targets;
    copy(mtype->begin(), mtype->end(), back_inserter(targets));

    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(targets.begin());
         itr != targets.end(); ++itr, ++idx)
    {
        const SpatiocyteWorld::coordinate_type source((*itr).coordinate);
        if (world_->get_voxel_pool_at(source) != mtype)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }
        const Integer rnd(rng->uniform_int(0, world_->num_neighbors(source)-1));
        const SpatiocyteWorld::coordinate_type destination(
                world_->get_neighbor_boundary(source, rnd));
        if (world_->can_move(source, destination))
        {
            if (rng->uniform(0,1) <= alpha)
                world_->move(source, destination, /*candidate=*/idx);
        }
        else
        {
            attempt_reaction_(*itr, destination, alpha);
        }
    }
}

void StepEvent::walk_on_surface_(const MoleculePool* mtype, const Real& alpha)
{
    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MoleculePool::container_type targets;
    copy(mtype->begin(), mtype->end(), back_inserter(targets));

    const VoxelPool* location(mtype->location());
    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(targets.begin());
         itr != targets.end(); ++itr, ++idx)
    {
        const SpatiocyteWorld::coordinate_type source((*itr).coordinate);
        if (world_->get_voxel_pool_at(source) != mtype)
        {
            // should skip if a voxel is not the target species.
            // when reaction has occured before, a voxel can be changed.
            continue;
        }

        const Integer num_neighbors(world_->num_neighbors(source));
        std::vector<unsigned int> nids;
        for (unsigned int i(0); i < num_neighbors; ++i)
            nids.push_back(i);
        ecell4::shuffle(*(rng.get()), nids);

        for (std::vector<unsigned int>::const_iterator nitr(nids.begin());
             nitr != nids.end(); ++nitr)
        {
            const SpatiocyteWorld::coordinate_type destination(
                    world_->get_neighbor_boundary(source, *nitr));
            const VoxelPool* target(world_->get_voxel_pool_at(destination));

            if (target->get_dimension() > mtype->get_dimension())
                continue;

            if (world_->can_move(source, destination))
            {
                if (rng->uniform(0,1) <= alpha)
                    world_->move(source, destination, /*candidate=*/idx);
            }
            else
            {
                attempt_reaction_(*itr, destination, alpha);
            }
            break;
        }
    }
}

std::pair<StepEvent::attempt_reaction_result_type, StepEvent::reaction_type>
StepEvent::attempt_reaction_(const SpatiocyteWorld::coordinate_id_pair_type& info,
                             const SpatiocyteWorld::coordinate_type to_coord,
                             const Real& alpha)
{
    const VoxelPool* src_vp(world_->get_voxel_pool_at(info.coordinate));
    const VoxelPool* dst_vp(world_->get_voxel_pool_at(to_coord));

    if (dst_vp->is_vacant())
    {
        return std::make_pair(NO_REACTION, reaction_type());
    }

    const Species& speciesA(src_vp->species());
    const Species& speciesB(dst_vp->species());

    const std::vector<ReactionRule> rules(
        model_->query_reaction_rules(speciesA, speciesB));

    if (rules.empty())
    {
        return std::make_pair(NO_REACTION, reaction_type());
    }

    const Real factor(calculate_dimensional_factor(src_vp, dst_vp,
                boost::const_pointer_cast<const SpatiocyteWorld>(world_)));

    const Real rnd(world_->rng()->uniform(0,1));
    Real accp(0.0);
    for (std::vector<ReactionRule>::const_iterator itr(rules.begin());
         itr != rules.end(); ++itr)
    {
        accp += (*itr).k() * factor * alpha;
        if (accp > 1)
        {
            std::cerr << "The total acceptance probability [" << accp
                << "] exceeds 1 for '" << speciesA.serial()
                << "' and '" << speciesB.serial() << "'." << std::endl;
        }
        if (accp >= rnd)
        {
            ReactionInfo rinfo(apply_second_order_reaction(
                        world_, *itr,
                        world_->make_pid_voxel_pair(src_vp, info),
                        world_->make_pid_voxel_pair(dst_vp, to_coord)));
            if (rinfo.has_occurred())
            {
                reaction_type reaction(std::make_pair(*itr, rinfo));
                push_reaction(reaction);
                return std::make_pair(REACTION_SUCCEEDED, reaction);
            }
            return std::make_pair(REACTION_FAILED, std::make_pair(*itr, rinfo));
        }
    }
    return std::make_pair(REACTION_FAILED, reaction_type());
}

} // spatiocyte

} // ecell4
