#include "SpatiocyteEvent.hpp"

namespace ecell4 {

namespace spatiocyte {

typedef SpatiocyteWorld::coordinate_type coord_type;

inline void
react_a2b(boost::shared_ptr<SpatiocyteWorld>& world,
          ReactionInfo& rinfo,
          const ReactionInfo::identified_voxel& voxelA,
          const Species& speciesB)
{
    const coord_type coordA(voxelA.second.coordinate());

    rinfo.add_reactant(voxelA);
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
       const ReactionInfo::identified_voxel& voxel)
{
    const coord_type coord(voxel.second.coordinate());

    rinfo.add_reactant(voxel);
    world->remove_voxel(coord);
    rinfo.add_product(world->get_voxel_at(coord));
}


/*
 * VanishmentEvent
 */

ReactionInfo VanishmentEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    ReactionInfo rinfo(world_->t());
    rinfo.add_reactant(voxel);
    world_->remove_voxel(voxel.second.coordinate());
    return rinfo;
}


/*
 * RearrangementEvent
 */

ReactionInfo RearrangementEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    const std::string serialLocationReactant(world_->get_molecule_info(reactant_).loc);
    const std::string serialLocationProduct(world_->get_molecule_info(product_).loc);

    ReactionInfo rinfo(world_->t());

    if (serialLocationReactant == serialLocationProduct)
    {
        react_a2b(world_, rinfo, voxel, product_);
    }
    else if (reactant_.serial() == serialLocationProduct)
    {
        /* XXX: deprecated
         *
         * should use ZerothOrderReactionEvent
         */
        rinfo.add_reactant(voxel);
        rinfo.add_product(world_->new_voxel(product_, voxel.second.coordinate()).first);
    }
    else if (serialLocationReactant == product_.serial())
    {
        /* XXX: deprecated
         *
         * should use VanishmentEvent
         */
        vanish(world_, rinfo, voxel);
    }
    else
    {
        const coord_type coord(voxel.second.coordinate());
        if (boost::optional<coord_type> neighbor
                = world_->check_neighbor(coord, serialLocationProduct))
        {
            rinfo.add_reactant(voxel);
            world_->remove_voxel(coord);
            rinfo.add_product(world_->new_voxel(product_, *neighbor).first);
        }
        else
        {
            return ReactionInfo(world_->t());
        }
    }

    push_product(product_);
    return rinfo;
}


/*
 * EliminationEvent
 */

ReactionInfo EliminationEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    // ReactionInfo rinfo(apply_a2bc(world_, voxel, product1_, product2_));

    const std::string serialReactant(reactant_.serial());
    const std::string serialLocationReactant(world_->get_molecule_info(reactant_).loc);
    const std::string serialLocationProduct1(world_->get_molecule_info(product1_).loc);
    const std::string serialLocationProduct2(world_->get_molecule_info(product2_).loc);

    ReactionInfo rinfo(world_->t());

    if (serialReactant == product1_.serial())
    {
        if (boost::optional<coord_type> neighbor
                = world_->check_neighbor(voxel.second.coordinate(), serialLocationProduct2))
        {
            rinfo.add_reactant(voxel);
            rinfo.add_product(voxel);
            rinfo.add_product(world_->new_voxel(product2_, *neighbor).first);
        }
    }
    else if (serialReactant == product2_.serial())
    {
        if (boost::optional<coord_type> neighbor
                = world_->check_neighbor(voxel.second.coordinate(), serialLocationProduct1))
        {
            rinfo.add_reactant(voxel);
            rinfo.add_product(voxel);
            rinfo.add_product(world_->new_voxel(product1_, *neighbor).first);
        }
    }
    else if (serialLocationReactant == serialLocationProduct1)
    {
        if (boost::optional<coord_type> neighbor
                = world_->check_neighbor(voxel.second.coordinate(), serialLocationProduct2))
        {
            react_a2b(world_, rinfo, voxel, product1_);
            rinfo.add_product(world_->new_voxel(product2_, *neighbor).first);
        }
    }
    else if (serialLocationReactant == serialLocationProduct2)
    {
        if (boost::optional<coord_type> neighbor
                = world_->check_neighbor(voxel.second.coordinate(), serialLocationProduct1))
        {
            react_a2b(world_, rinfo, voxel, product2_);
            rinfo.add_product(world_->new_voxel(product1_, *neighbor).first);
        }
    }
    else if (serialLocationReactant == product1_.serial())
    {
        /* XXX: deprecated
         *
         * should use RearrangementEvent
         */
        if (boost::optional<coord_type> neighbor
                = world_->check_neighbor(voxel.second.coordinate(), serialLocationProduct2))
        {
            vanish(world_, rinfo, voxel);
            rinfo.add_product(world_->new_voxel(product2_, *neighbor).first);
        }
    }
    else if (serialLocationReactant == product2_.serial())
    {
        /* XXX: deprecated
         *
         * should use RearrangementEvent
         */
        if (boost::optional<coord_type> neighbor
                = world_->check_neighbor(voxel.second.coordinate(), serialLocationProduct1))
        {
            rinfo.add_product(world_->new_voxel(product1_, *neighbor).first);
            vanish(world_, rinfo, voxel);
        }
    }
    else if (serialReactant == serialLocationProduct1)
    {
        const coord_type coord(voxel.second.coordinate());
        if (boost::optional<coord_type> neighbor
                = world_->check_neighbor(coord, serialLocationProduct2))
        {
            rinfo.add_reactant(voxel);
            rinfo.add_product(world_->new_voxel(product1_, coord).first);
            rinfo.add_product(world_->new_voxel(product2_, *neighbor).first);
        }
    }
    else if (serialReactant == serialLocationProduct2)
    {
        const coord_type coord(voxel.second.coordinate());
        if (boost::optional<coord_type> neighbor
                = world_->check_neighbor(coord, serialLocationProduct1))
        {
            rinfo.add_reactant(voxel);
            rinfo.add_product(world_->new_voxel(product1_, *neighbor).first);
            rinfo.add_product(world_->new_voxel(product2_, coord).first);
        }
    }

    push_product(product1_);
    push_product(product2_);

    return rinfo;
}

} // spatiocyte

} // ecell4
