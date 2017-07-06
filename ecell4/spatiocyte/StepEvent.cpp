#include "SpatiocyteEvent.hpp"
#include "utils.hpp"

namespace ecell4
{

namespace spatiocyte
{

inline
ReactionInfo
apply_vanishment(boost::shared_ptr<SpatiocyteWorld> world,
                 const ReactionInfo::identified_voxel& p0,
                 const ReactionInfo::identified_voxel& p1)
{
    world->remove_voxel(p0.second.coordinate());
    world->remove_voxel(p1.second.coordinate());

    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(p0);
    rinfo.add_reactant(p1);

    return rinfo;
}

inline
ReactionInfo
apply_ab2c(boost::shared_ptr<SpatiocyteWorld> world,
           const ReactionInfo::identified_voxel& reactantA,
           const ReactionInfo::identified_voxel& reactantB,
           const Species& product_species)
{
    const std::string& serialA(get_serial(world, reactantA.second.coordinate()));
    const std::string& serialB(get_serial(world, reactantB.second.coordinate()));

    const std::string& locationA(get_location(world, reactantA.second.coordinate()));
    const std::string& locationB(get_location(world, reactantB.second.coordinate()));
    const std::string& locationC(world->get_molecule_info(product_species).loc);

    SpatiocyteWorld::coordinate_type new_coordinate;

    if (locationB == locationC)
    {
        world->remove_voxel(reactantB.second.coordinate());
        world->remove_voxel(reactantA.second.coordinate());

        new_coordinate = reactantB.second.coordinate();
    }
    else if (locationA == locationC)
    {
        world->remove_voxel(reactantA.second.coordinate());
        world->remove_voxel(reactantB.second.coordinate());

        new_coordinate = reactantA.second.coordinate();
    }
    else if (serialB == locationC)
    {
        world->remove_voxel(reactantA.second.coordinate());

        new_coordinate = reactantB.second.coordinate();
    }
    else if (serialA == locationC)
    {
        world->remove_voxel(reactantB.second.coordinate());

        new_coordinate = reactantA.second.coordinate();
    }
    else
    {
        return ReactionInfo(world->t());
    }

    std::pair<ReactionInfo::identified_voxel, bool> new_mol(world->new_voxel(product_species, new_coordinate));

    ReactionInfo rinfo(world->t());

    rinfo.add_reactant(reactantA);
    rinfo.add_reactant(reactantB);
    rinfo.add_product(new_mol.first);

    return rinfo;
}

inline
ReactionInfo
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
            return apply_ab2c(world, p0, p1, products.at(0));
        case 2:
            return apply_ab2cd(world, p0, p1, products.at(0), products.at(1));
        default:
            return ReactionInfo(world->t());
    }
}

StepEvent::StepEvent(boost::shared_ptr<Model> model,
                     boost::shared_ptr<SpatiocyteWorld> world,
                     const Species& species, const Real& t, const Real alpha)
    : SpatiocyteEvent(t), model_(model), world_(world), species_(species), alpha_(alpha)
{
    const MoleculeInfo minfo(world_->get_molecule_info(species));
    const Real R(world_->voxel_radius());
    const Real D(minfo.D);
    const Real sqRperD(pow(R, 2.0)/D);
    const VoxelPool* vpool(world_->find_voxel_pool(species));
    if (D <= 0)
    {
        dt_ = inf;
    } else if(vpool->get_dimension() == Shape::THREE) {
        dt_ = 2 * sqRperD / 3 * alpha_;
    } else if(vpool->get_dimension() == Shape::TWO) {
        // TODO: Regular Lattice
        // dt_  = pow((2*sqrt(2.0)+4*sqrt(3.0)+3*sqrt(6.0)+sqrt(22.0))/
        //           (6*sqrt(2.0)+4*sqrt(3.0)+3*sqrt(6.0)), 2) * sqRperD * alpha_;
        dt_ = sqRperD * alpha_;
    } else if(vpool->get_dimension() == Shape::ONE) {
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

    const MoleculePool* mpool(world_->find_molecule_pool(species_));

    if (mpool->get_dimension() == Shape::THREE)
        walk_in_space_(mpool, alpha);
    else // dimension == TWO, etc.
        walk_on_surface_(mpool, alpha);
}

void StepEvent::walk_in_space_(const MoleculePool* mpool, const Real& alpha)
{
    const boost::shared_ptr<RandomNumberGenerator> rng(world_->rng());
    MoleculePool::container_type targets;
    copy(mpool->begin(), mpool->end(), back_inserter(targets));

    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(targets.begin());
         itr != targets.end(); ++itr, ++idx)
    {
        const SpatiocyteWorld::coordinate_type source((*itr).coordinate);

        // skip when the voxel is not the target species.
        // former reactions may change the voxel.
        if (world_->get_voxel_pool_at(source) != mpool)
            continue;

        const Integer rnd(rng->uniform_int(0, world_->num_neighbors(source)-1));
        const SpatiocyteWorld::coordinate_type destination(world_->get_neighbor_boundary(source, rnd));

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

void StepEvent::walk_on_surface_(const MoleculePool* mpool, const Real& alpha)
{
    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MoleculePool::container_type targets;
    copy(mpool->begin(), mpool->end(), back_inserter(targets));

    const VoxelPool* location(mpool->location());
    const Shape::dimension_kind dimension(mpool->get_dimension());

    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(targets.begin());
         itr != targets.end(); ++itr, ++idx)
    {
        const SpatiocyteWorld::coordinate_type source((*itr).coordinate);

        // skip when the voxel is not the target species.
        // former reactions may change the voxel.
        if (world_->get_voxel_pool_at(source) != mpool)
            continue;

        // list up the neighbors whose location is ok.
        std::vector<SpatiocyteWorld::coordinate_type> neighbors;
        for (unsigned int i(0); i < world_->num_neighbors(source); ++i)
        {
            const SpatiocyteWorld::coordinate_type
                neighbor(world_->get_neighbor_boundary(source, i));

            if (world_->get_voxel_pool_at(neighbor)->get_dimension() <= dimension)
                neighbors.push_back(neighbor);
        }

        if (neighbors.size() == 0)
            continue;

        const SpatiocyteWorld::coordinate_type
            destination(neighbors.at(rng->uniform_int(0, neighbors.size()-1)));

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

std::pair<StepEvent::attempt_reaction_result_type, StepEvent::reaction_type>
StepEvent::attempt_reaction_(const SpatiocyteWorld::coordinate_id_pair_type& info,
                             SpatiocyteWorld::coordinate_type to_coord,
                             const Real& alpha)
{
    to_coord = world_->get_adjoining_or_self(to_coord);

    const VoxelPool* vpA(world_->get_voxel_pool_at(info.coordinate));
    const VoxelPool* vpB(world_->get_voxel_pool_at(to_coord));

    if (vpB->is_vacant())
    {
        return std::make_pair(NO_REACTION, reaction_type());
    }

    const Species& speciesA(vpA->species());
    const Species& speciesB(vpB->species());

    const std::vector<ReactionRule> rules(model_->query_reaction_rules(speciesA, speciesB));

    if (rules.empty())
    {
        return std::make_pair(NO_REACTION, reaction_type());
    }

    const Real factor(calculate_dimensional_factor(vpA, vpB,
                boost::const_pointer_cast<const SpatiocyteWorld>(world_)));

    const Real rnd(world_->rng()->uniform(0,1));
    Real accp(0.0);
    for (std::vector<ReactionRule>::const_iterator rritr(rules.begin());
         rritr != rules.end(); ++rritr)
    {
        const ReactionRule &reaction_rule(*rritr);
        accp += reaction_rule.k() * factor * alpha;
        if (accp > 1)
        {
            std::cerr << "The total acceptance probability [" << accp
                << "] exceeds 1 for '" << speciesA.serial()
                << "' and '" << speciesB.serial() << "'." << std::endl;
        }
        if (accp >= rnd)
        {
            ReactionInfo rinfo(apply_second_order_reaction(
                        world_, reaction_rule,
                        world_->make_pid_voxel_pair(vpA, info),
                        world_->make_pid_voxel_pair(vpB, to_coord)));
            if (rinfo.has_occurred())
            {
                reaction_type reaction(std::make_pair(reaction_rule, rinfo));
                push_reaction(reaction);
                return std::make_pair(REACTION_SUCCEEDED, reaction);
            }
            return std::make_pair(REACTION_FAILED, std::make_pair(reaction_rule, rinfo));
        }
    }
    return std::make_pair(REACTION_FAILED, reaction_type());
}

} // spatiocyte

} // ecell4
