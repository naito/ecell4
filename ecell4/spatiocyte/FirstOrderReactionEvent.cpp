#include "SpatiocyteEvent.hpp"

namespace ecell4 {

namespace spatiocyte {

inline ReactionInfo
apply_a2b(boost::shared_ptr<SpatiocyteWorld> world,
          const ReactionInfo::identified_voxel& p,
          const Species& product_species)
{
    const SpatiocyteWorld::coordinate_type coord(p.second.coordinate());
    const std::string bloc(world->get_molecule_info(product_species).loc);
    const std::string aserial(get_serial(world, coord));
    const std::string aloc(get_location(world, coord));
    const std::string bserial(product_species.serial());

    ReactionInfo rinfo(world->t());

    if (aserial == bloc || aloc == bloc || aloc == bserial)
    {
        // A is the location of B (B can be placed on A),
        // or A is on the location of B,
        // or A is on B.
        rinfo.add_reactant(p);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of B
            world->remove_voxel(p.second.coordinate());
        }

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world->new_voxel(product_species, coord));
            rinfo.add_product(new_mol.first);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            rinfo.add_product(world->get_voxel_at(coord));
        }
    }
    else
    {
        // A is NOT on the location of B.
        // B must be released into a neighbor, which is the location of B
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(coord, bloc));

        if (neighbor.second)
        {
            // The neighbor is the location of B.
            // Place B at the neighbor, and remove A.
            rinfo.add_reactant(p);

            world->remove_voxel(p.second.coordinate());
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol(
                world->new_voxel(product_species, neighbor.first));

            rinfo.add_product(new_mol.first);
        }
    }
    return rinfo;
}

inline ReactionInfo
apply_a2bc(boost::shared_ptr<SpatiocyteWorld> world,
           const ReactionInfo::identified_voxel& p,
           const Species& product_species0,
           const Species& product_species1)
{
    // A (pinfo) becomes B and C (product_species0 and product_species1)
    // At least, one of A and B must be placed at the neighbor.
    const SpatiocyteWorld::coordinate_type coord(p.second.coordinate());
    const std::string
        bserial(product_species0.serial()),
        cserial(product_species1.serial()),
        bloc(world->get_molecule_info(product_species0).loc),
        cloc(world->get_molecule_info(product_species1).loc);
    const std::string aserial(get_serial(world, coord));
    const std::string aloc(get_location(world, coord));

    ReactionInfo rinfo(world->t());

    if (aserial == bloc || aloc == bloc || aloc == bserial)
    {
        // C must be placed at the neighbor

        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(coord, cloc));
        const std::string nserial(get_serial(world, neighbor.first));
        const std::string nloc(get_location(world, neighbor.first));

        if (!neighbor.second)
        {
            //TODO: C cannot be on the neighbor.
            return rinfo;
        }

        rinfo.add_reactant(p);

        if (aserial != bloc)
        {
            // Remove A once if A is not the location of a new B-molecule
            world->remove_voxel(p.second.coordinate());
        }

        // No need to remove the neighbor because it's the location of C
        // world->remove_voxel(neighbor.first);

        if (aloc != bserial)
        {
            // Place a new B-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
                    world->new_voxel(product_species0, coord));
            rinfo.add_product(new_mol0.first);
        }
        else
        {
            // When B is the location of A, it's enough to remove A
            rinfo.add_product(world->get_voxel_at(coord));
        }

        // Place a new C-molecule at the neighbor
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
            world->new_voxel(product_species1, neighbor.first));
        rinfo.add_product(new_mol1.first);
        return rinfo;
    }
    else if (aserial == cloc || aloc == cloc || aloc == cserial)
    {
        // A is the locaiton of C,
        // or A is on the location of C,
        // or C is the location of A
        // B must be placed at the neighbor
        std::pair<SpatiocyteWorld::coordinate_type, bool>
            neighbor(world->check_neighbor(coord, bloc));
        const std::string nserial(get_serial(world, neighbor.first));
        const std::string nloc(get_location(world, neighbor.first));

        if (!neighbor.second)
        {
            //TODO: B cannot be on the neighbor.
            return rinfo;
        }

        rinfo.add_reactant(p);

        if (aserial != cloc)
        {
            // Remove A once if A is not the location of a new C-molecule
            world->remove_voxel(p.second.coordinate());
        }

        // No need to remove the neighbor because it's the location of B
        // world->remove_voxel(neighbor.first);

        // Place a new B-molecule at the neighbor
        std::pair<std::pair<ParticleID, Voxel>, bool> new_mol0(
            world->new_voxel(product_species0, neighbor.first));
        rinfo.add_product(new_mol0.first);

        if (aloc != cserial)
        {
            // Place a new C-molecule at the position of A
            std::pair<std::pair<ParticleID, Voxel>, bool> new_mol1(
                world->new_voxel(product_species1, coord));
            rinfo.add_product(new_mol1.first);
        }
        else
        {
            // When C is the location of A, it's enough to remove A
            rinfo.add_product(world->get_voxel_at(coord));
        }
        return rinfo;
    }
    return rinfo;
}

FirstOrderReactionEvent::FirstOrderReactionEvent(
    boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(t), world_(world), rule_(rule)
{
    //assert(rule_.reactants().size() == 1);
    time_ = t + draw_dt();
}

void FirstOrderReactionEvent::fire_()
{
    const ReactionInfo::identified_voxel& p(
            world_->choice(*(rule_.reactants().begin())));
    const ReactionRule::product_container_type& products(rule_.products());

    time_ += draw_dt();
    switch (products.size())
    {
        case 0:
            {
                world_->remove_voxel(p.second.coordinate());
                ReactionInfo rinfo(world_->t());
                rinfo.add_reactant(p);
                push_reaction(std::make_pair(rule_, rinfo));
            }
            break;
        case 1:
            push_reaction(std::make_pair(rule_, apply_a2b(world_, p, *(products.begin()))));
            break;
        case 2:
            {
                ReactionInfo rinfo(apply_a2bc(world_, p,
                            *(products.begin()), (*(++products.begin()))));
                if (rinfo.has_occurred())
                    push_reaction(std::make_pair(rule_, rinfo));
            }
            break;
    }
}

Real FirstOrderReactionEvent::draw_dt()
{
    const Species& reactant(*(rule_.reactants().begin()));
    const Integer num_r(world_->num_voxels_exact(reactant));
    const Real k(rule_.k());
    const Real p = k * num_r;
    Real dt(inf);
    if (p != 0.)
    {
        const Real rnd(world_->rng()->uniform(0.,1.));
        dt = - log(1 - rnd) / p;
    }
    return dt;
}

} // spatiocyte

} // ecell4
