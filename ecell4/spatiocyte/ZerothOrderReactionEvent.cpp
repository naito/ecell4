#include "SpatiocyteEvent.hpp"

namespace ecell4 {

namespace spatiocyte {

typedef SpatiocyteWorld::coordinate_type coord_type;

ZerothOrderReactionEvent::ZerothOrderReactionEvent(
    boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t)
    : SpatiocyteEvent(world, t), rule_(rule)
{
    time_ = t + draw_dt();
}

void ZerothOrderReactionEvent::fire_()
{
    ReactionInfo rinfo(world_->t());

    std::vector<coord_type> coordinates(world_->inner_size());
    for (coord_type i(0); i < world_->inner_size(); ++i)
        coordinates.push_back(i);

    for (ReactionRule::product_container_type::const_iterator i(rule_.products().begin());
        i != rule_.products().end(); ++i)
    {
        const Species& sp(*i);
        const MoleculeInfo& info(world_->get_molecule_info(sp));

        push_product(sp);

        shuffle(*(world_->rng()), coordinates);
        for (std::vector<coord_type>::const_iterator itr(coordinates.begin());
             itr != coordinates.end(); ++itr)
        {
            const coord_type coord(world_->inner2coordinate(*itr));
            const Voxel v(sp, coord, info.radius, info.D, info.loc);

            if (world_->on_structure(v))
                continue;

            const std::pair<std::pair<ParticleID, Voxel>, bool> retval(world_->new_voxel(v));
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
    if (p > 0)
    {
        const Real rnd(world_->rng()->uniform(0.,1.));
        dt = - log(1 - rnd) / p;
    }
    return dt;
}

} // spatiocyte

} // ecell4
