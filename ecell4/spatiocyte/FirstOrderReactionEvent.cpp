#include "SpatiocyteEvent.hpp"

namespace ecell4 {

namespace spatiocyte {

typedef SpatiocyteWorld::coordinate_type coord_type;

inline void
react_a2b(boost::shared_ptr<SpatiocyteWorld>& world,
          ReactionInfo& rinfo,
          const coord_type coordA,
          const Species& speciesB)
{
    world->remove_voxel(coordA);
    rinfo.add_product(world->new_voxel(speciesB, coordA).first);
}

/* Deprecated function
 *
 * If this function is needed, the reaction is invalid.
 */
inline void
vanish(boost::shared_ptr<SpatiocyteWorld>& world,
       ReactionInfo& rinfo,
       const coord_type coord)
{
    world->remove_voxel(coord);
    rinfo.add_product(world->get_voxel_at(coord));
}

ReactionInfo
apply_a2b(boost::shared_ptr<SpatiocyteWorld> world,
          const ReactionInfo::identified_voxel& voxelA,
          const Species& speciesB)
{
    const coord_type coordA(voxelA.second.coordinate());
    const std::string serialA(get_serial(world, coordA)),
                      serialB(speciesB.serial()),
                      serialLocationA(get_location(world, coordA)),
                      serialLocationB(world->get_molecule_info(speciesB).loc);

    ReactionInfo rinfo(world->t());
    rinfo.add_reactant(voxelA);

    if (serialLocationA == serialLocationB)
        react_a2b(world, rinfo, coordA, speciesB);
    else if (serialA == serialLocationB)
    {
        /* XXX: deprecated
         *
         * have to use the 0th order reaction not the 1st order reaction
         */
        rinfo.add_product(world->new_voxel(speciesB, coordA).first);
    }
    else if (serialLocationA == serialB)
    {
        /* XXX: deprecated
         *
         * have to use the vanishment reaction not the a2b reaction.
         */
        vanish(world, rinfo, coordA);
    }
    else
    {
        std::pair<coord_type, bool> neighbor(world->check_neighbor(coordA, serialLocationB));

        if (neighbor.second)
        {
            rinfo.add_reactant(voxelA);
            world->remove_voxel(coordA);
            rinfo.add_product(world->new_voxel(speciesB, neighbor.first).first);
        }
        else
            return ReactionInfo(world->t());
    }
    return rinfo;
}

ReactionInfo
apply_a2bc(boost::shared_ptr<SpatiocyteWorld> world,
           const ReactionInfo::identified_voxel& voxelA,
           const Species& speciesB, const Species& speciesC)
{
    const coord_type coordA(voxelA.second.coordinate());
    const std::string serialA(get_serial(world, coordA)),
                      serialB(speciesB.serial()),
                      serialC(speciesC.serial()),
                      serialLocationA(get_location(world, coordA)),
                      serialLocationB(world->get_molecule_info(speciesB).loc),
                      serialLocationC(world->get_molecule_info(speciesC).loc);

    ReactionInfo rinfo(world->t());

    if (serialA == serialB)
    {
        std::pair<coord_type, bool> neighbor(world->check_neighbor(coordA, serialLocationC));
        if (neighbor.second)
        {
            rinfo.add_reactant(voxelA);
            rinfo.add_product(voxelA);
            rinfo.add_product(world->new_voxel(speciesC, neighbor.first).first);
        }
    }
    else if (serialA == serialC)
    {
        std::pair<coord_type, bool> neighbor(world->check_neighbor(coordA, serialLocationB));
        if (neighbor.second)
        {
            rinfo.add_reactant(voxelA);
            rinfo.add_product(voxelA);
            rinfo.add_product(world->new_voxel(speciesB, neighbor.first).first);
        }
    }
    else if (serialLocationA == serialLocationB)
    {
        std::pair<coord_type, bool> neighbor(world->check_neighbor(coordA, serialLocationC));
        if (neighbor.second)
        {
            rinfo.add_reactant(voxelA);
            react_a2b(world, rinfo, coordA, speciesB);
            rinfo.add_product(world->new_voxel(speciesC, neighbor.first).first);
        }
    }
    else if (serialLocationA == serialLocationC)
    {
        std::pair<coord_type, bool> neighbor(world->check_neighbor(coordA, serialLocationB));
        if (neighbor.second)
        {
            rinfo.add_reactant(voxelA);
            rinfo.add_product(world->new_voxel(speciesB, neighbor.first).first);
            react_a2b(world, rinfo, coordA, speciesC);
        }
    }
    else if (serialLocationA == serialB)
    {
        /* XXX: deprecated
         *
         * have to use the a2b reaction not the a2bc reaction.
         */
        std::pair<coord_type, bool> neighbor(world->check_neighbor(coordA, serialLocationC));
        if (neighbor.second)
        {
            rinfo.add_reactant(voxelA);
            vanish(world, rinfo, coordA);
            rinfo.add_product(world->new_voxel(speciesC, neighbor.first).first);
        }
    }
    else if (serialLocationA == serialC)
    {
        /* XXX: deprecated
         *
         * have to use the a2b reaction not the a2bc reaction.
         */
        std::pair<coord_type, bool> neighbor(world->check_neighbor(coordA, serialLocationB));
        if (neighbor.second)
        {
            rinfo.add_reactant(voxelA);
            rinfo.add_product(world->new_voxel(speciesB, neighbor.first).first);
            vanish(world, rinfo, coordA);
        }
    }
    else if (serialA == serialLocationB)
    {
        /* XXX: deprecated
         *
         * B and C arise simultaneously. This is invalid reaction.
         */
        std::pair<coord_type, bool> neighbor(world->check_neighbor(coordA, serialLocationC));
        if (neighbor.second)
        {
            rinfo.add_reactant(voxelA);
            rinfo.add_product(world->new_voxel(speciesB, coordA).first);
            rinfo.add_product(world->new_voxel(speciesC, neighbor.first).first);
        }
    }
    else if (serialA == serialLocationC)
    {
        /* XXX: deprecated
         *
         * B and C arise simultaneously. This is invalid reaction.
         */
        std::pair<coord_type, bool> neighbor(world->check_neighbor(coordA, serialLocationB));
        if (neighbor.second)
        {
            rinfo.add_reactant(voxelA);
            rinfo.add_product(world->new_voxel(speciesB, neighbor.first).first);
            rinfo.add_product(world->new_voxel(speciesC, coordA).first);
        }
    }

    return rinfo;
}

FirstOrderReactionEvent::FirstOrderReactionEvent(boost::shared_ptr<SpatiocyteWorld> world,
                                                 const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(world, t), rule_(rule)
{
    //assert(rule_.reactants().size() == 1);
    time_ = t + draw_dt();
}

void FirstOrderReactionEvent::fire_()
{
    if (world_->num_voxels_exact(reactant_()) == 0)
    {
        time_ += draw_dt();
        return;
    }

    const ReactionInfo::identified_voxel& voxel(world_->choice(reactant_()));
    const ReactionRule::product_container_type& products(rule_.products());

    switch (products.size())
    {
        case 0:
            {
                ReactionInfo rinfo(world_->t());
                rinfo.add_reactant(voxel);
                world_->remove_voxel(voxel.second.coordinate());
                push_reaction(std::make_pair(rule_, rinfo));
            }
            break;
        case 1:
            {
                push_product(products.at(0));
                ReactionInfo rinfo(apply_a2b(world_, voxel, products.at(0)));
                if (rinfo.has_occurred())
                    push_reaction(std::make_pair(rule_, rinfo));
            }
            break;
        case 2:
            {
                push_product(products.at(0));
                push_product(products.at(1));
                ReactionInfo rinfo(apply_a2bc(world_, voxel, products.at(0), products.at(1)));
                if (rinfo.has_occurred())
                {
                    push_reaction(std::make_pair(rule_, rinfo));
                }
            }
            break;
    }
    time_ += draw_dt();
}

Real FirstOrderReactionEvent::draw_dt()
{
    const Integer num(world_->num_voxels_exact(reactant_()));
    const Real k(rule_.k());
    const Real p = k * num;
    Real dt(inf);
    if (p != 0.)
    {
        const Real rnd(world_->rng()->uniform(0.,1.));
        dt = - log(1 - rnd) / p;
    }
    prev_num_voxels_ = num;
    return dt;
}

} // spatiocyte

} // ecell4
