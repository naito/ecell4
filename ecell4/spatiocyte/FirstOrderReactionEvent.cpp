#include "FirstOrderReactionEvent.hpp"

namespace ecell4 {

namespace spatiocyte {

typedef SpatiocyteWorld::coordinate_type coord_type;

static inline FirstOrderReactionEvent*
generate_a2b_reaction_event(boost::shared_ptr<SpatiocyteWorld> world,
                            const ReactionRule& rule,
                            const Real& t)
{
    const Species& reactant(rule.reactants().at(0));
    const Species& product(rule.products().at(0));

    const std::string serialLocationReactant(world->get_molecule_info(reactant).loc);
    const std::string serialLocationProduct(world->get_molecule_info(product).loc);

    if (serialLocationReactant == serialLocationProduct)
        return new RearrangementEvent(world, rule, t);

    if (reactant.serial() == serialLocationProduct)
        return new GenerationEvent(world, rule, t);

    if (serialLocationReactant == product.serial())
        return new VanishmentEvent(world, rule, t);

    return new DesorptionEvent(world, rule, product, serialLocationProduct, t);

}

static inline FirstOrderReactionEvent*
generate_a2bc_reaction_event(boost::shared_ptr<SpatiocyteWorld> world,
                             const ReactionRule& rule,
                             const Real& t)
{
    const Species& reactant(rule.reactants().at(0));
    const Species& product1(rule.products().at(0));
    const Species& product2(rule.products().at(1));

    const std::string serialLocationReactant(world->get_molecule_info(reactant).loc);
    const std::string serialLocationProduct1(world->get_molecule_info(product1).loc);
    const std::string serialLocationProduct2(world->get_molecule_info(product2).loc);

    if (reactant.serial() == product1.serial())
        return new GenerationBesideEvent(world, rule, product2, serialLocationProduct2, t);

    if (reactant.serial() == product2.serial())
        return new GenerationBesideEvent(world, rule, product1, serialLocationProduct1, t);

    if (serialLocationReactant == serialLocationProduct1)
        return new EliminationEvent(world, rule, product1, product2, serialLocationProduct2, t);

    if (serialLocationReactant == serialLocationProduct2)
        return new EliminationEvent(world, rule, product2, product1, serialLocationProduct1, t);

    if (reactant.serial() == serialLocationProduct1)
        return new GenerationTwoEvent(world, rule, product1, product2, serialLocationProduct2, t);

    if (reactant.serial() == serialLocationProduct2)
        return new GenerationTwoEvent(world, rule, product2, product1, serialLocationProduct1, t);

    if (serialLocationReactant == product1.serial())
        return new DesorptionEvent(world, rule, product2, serialLocationProduct2, t);

    if (serialLocationReactant == product2.serial())
        return new DesorptionEvent(world, rule, product1, serialLocationProduct1, t);

    throw NotSupported("The type of a given ReactionRule is not supported.");
}

FirstOrderReactionEvent*
generate_first_order_reaction_event(boost::shared_ptr<SpatiocyteWorld> world,
                                    const ReactionRule& rule,
                                    const Real& t)
{
    switch (rule.products().size())
    {
        case 0:
            return new VanishmentEvent(world, rule, t);
        case 1:
            return generate_a2b_reaction_event(world, rule, t);
        case 2:
            return generate_a2bc_reaction_event(world, rule, t);
    }
    throw NotSupported("The size of products is supported only for 0, 1 or 2.");
}

static inline void
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
static inline void
vanish(boost::shared_ptr<SpatiocyteWorld>& world,
       ReactionInfo& rinfo,
       const ReactionInfo::identified_voxel& voxel)
{
    const coord_type coord(voxel.second.coordinate());

    rinfo.add_reactant(voxel);
    world->remove_voxel(coord);
    rinfo.add_product(world->get_voxel_at(coord));
}


boost::optional<ReactionInfo>
VanishmentEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    ReactionInfo rinfo(world_->t());
    rinfo.add_reactant(voxel);
    world_->remove_voxel(voxel.second.coordinate());
    return rinfo;
}

boost::optional<ReactionInfo>
RearrangementEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    ReactionInfo rinfo(world_->t());
    react_a2b(world_, rinfo, voxel, product_);
    return rinfo;
}

boost::optional<ReactionInfo>
GenerationEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    ReactionInfo rinfo(world_->t());
    rinfo.add_reactant(voxel);
    rinfo.add_product(world_->new_voxel(product_, voxel.second.coordinate()).first);
    return rinfo;
}

boost::optional<ReactionInfo>
DesorptionEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    const coord_type coord(voxel.second.coordinate());
    if (boost::optional<coord_type> neighbor
            = world_->check_neighbor(coord, serial_location_product_))
    {
        ReactionInfo rinfo(world_->t());
        rinfo.add_reactant(voxel);
        world_->remove_voxel(coord);
        rinfo.add_product(world_->new_voxel(product_, *neighbor).first);
        return rinfo;
    }

    return boost::none;
}

boost::optional<ReactionInfo>
EliminationEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    if (boost::optional<coord_type> neighbor
            = world_->check_neighbor(voxel.second.coordinate(), serial_location_product_))
    {
        ReactionInfo rinfo(world_->t());
        react_a2b(world_, rinfo, voxel, product1_);
        rinfo.add_product(world_->new_voxel(product2_, *neighbor).first);
        return rinfo;
    }

    return boost::none;
}

boost::optional<ReactionInfo>
GenerationBesideEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    if (boost::optional<coord_type> neighbor
            = world_->check_neighbor(voxel.second.coordinate(), serial_location_product_))
    {
        ReactionInfo rinfo(world_->t());
        rinfo.add_reactant(voxel);
        rinfo.add_product(voxel);
        rinfo.add_product(world_->new_voxel(product_, *neighbor).first);
        return rinfo;
    }

    return boost::none;
}

boost::optional<ReactionInfo>
GenerationTwoEvent::react(const ReactionInfo::identified_voxel& voxel)
{
    const coord_type coord(voxel.second.coordinate());
    if (boost::optional<coord_type> neighbor
            = world_->check_neighbor(coord, serial_location_product_))
    {
        ReactionInfo rinfo(world_->t());
        rinfo.add_reactant(voxel);
        rinfo.add_product(world_->new_voxel(product1_, coord).first);
        rinfo.add_product(world_->new_voxel(product2_, *neighbor).first);
        return rinfo;
    }

    return boost::none;
}

} // spatiocyte

} // ecell4
