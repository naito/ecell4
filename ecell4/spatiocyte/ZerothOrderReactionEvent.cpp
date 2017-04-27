#include "SpatiocyteEvent.hpp"

namespace ecell4 {

namespace spatiocyte {

ZerothOrderReactionEvent::ZerothOrderReactionEvent(
    boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(t), world_(world), rule_(rule)
{
    time_ = t + draw_dt();
}

void ZerothOrderReactionEvent::fire_()
{
    ReactionInfo rinfo(world_->t());

    for (ReactionRule::product_container_type::const_iterator
        i(rule_.products().begin());
        i != rule_.products().end(); ++i)
    {
        const Species& sp(*i);
        const MoleculeInfo info(world_->get_molecule_info(sp));

        while (true) //TODO: Avoid an inifinite loop
        {
            // const SpatiocyteWorld::coordinate_type
            //     coord(world_->rng()->uniform_int(0, world_->size() - 1));
            const SpatiocyteWorld::coordinate_type
                coord(world_->inner2coordinate(
                            world_->rng()->uniform_int(0, world_->inner_size() - 1)));
            const Voxel v(sp, coord, info.radius, info.D, info.loc);

            if (world_->on_structure(v))
            {
                continue;
            }

            const std::pair<std::pair<ParticleID, Voxel>, bool>
                retval(world_->new_voxel(v));
            if (retval.second)
            {
                rinfo.add_product(retval.first);
                break;
            }
        }
    }
    push_reaction(std::make_pair(rule_, rinfo));
    time_ += draw_dt();
}

Real ZerothOrderReactionEvent::draw_dt()
{
    const Real k(rule_.k());
    const Real p = k * world_->volume();
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
