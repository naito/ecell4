#ifndef ECELL4_SPATIOCYTE_FIRST_ORDER_REACTION_EVENT_HPP
#define ECELL4_SPATIOCYTE_FIRST_ORDER_REACTION_EVENT_HPP

#include <boost/optional.hpp>
#include "SpatiocyteEvent.hpp"

namespace ecell4
{

namespace spatiocyte
{

struct FirstOrderReactionEvent : SpatiocyteEvent
{
    FirstOrderReactionEvent(boost::shared_ptr<SpatiocyteWorld> world,
                            const ReactionRule& rule,
                            const Real& t)
        : SpatiocyteEvent(world, t), rule_(rule), reactant_(rule_.reactants().at(0))
    {
        //assert(rule_.reactants().size() == 1);
        time_ = t + draw_dt();
    }

    virtual ~FirstOrderReactionEvent() {}

    virtual boost::optional<ReactionInfo> react(const ReactionInfo::identified_voxel& voxel) = 0;

    virtual void fire_()
    {
        if (world_->num_voxels_exact(reactant_) != 0)
        {
            if (boost::optional<ReactionInfo> rinfo = react(world_->choice(reactant_)))
            {
                push_reaction(std::make_pair(rule_, *rinfo));

                for (ReactionInfo::container_type::const_iterator itr(rinfo->products().begin());
                        itr != rinfo->products().end(); ++itr)
                {
                    push_product((*itr).second.species());
                }
            }
        }
        time_ += draw_dt();
    }

    Real draw_dt()
    {
        const Integer num(world_->num_voxels_exact(reactant_));
        const Real p = rule_.k() * num;
        Real dt(inf);
        if (p > 0.)
        {
            const Real rnd(world_->rng()->uniform(0.,1.));
            dt = - log(1 - rnd) / p;
        }
        prev_num_voxels_ = num;
        return dt;
    }

    virtual void interrupt(Real const& t)
    {
        if (world_->num_voxels_exact(reactant_) != prev_num_voxels_) {
            time_ = t + draw_dt();
        }
    }

protected:

    const ReactionRule rule_;
    const Species& reactant_;

private:

    Integer prev_num_voxels_;

};

struct VanishmentEvent : FirstOrderReactionEvent
{
    VanishmentEvent(boost::shared_ptr<SpatiocyteWorld> world,
                    const ReactionRule& rule,
                    const Real& t)
        : FirstOrderReactionEvent(world, rule, t)
    {
        // assert(rule_.products().size() == 0);
    }

    virtual boost::optional<ReactionInfo> react(const ReactionInfo::identified_voxel& voxel);
};


struct RearrangementEvent : FirstOrderReactionEvent
{
    RearrangementEvent(boost::shared_ptr<SpatiocyteWorld> world,
                       const ReactionRule& rule,
                       const Real& t)
        : FirstOrderReactionEvent(world, rule, t),
          product_(rule_.products().at(0))
    {
        // assert(rule_.products().size() == 1);
    }

    virtual boost::optional<ReactionInfo> react(const ReactionInfo::identified_voxel& voxel);

protected:

    const Species product_;

};

struct GenerationEvent : FirstOrderReactionEvent
{
    GenerationEvent(boost::shared_ptr<SpatiocyteWorld> world,
                    const ReactionRule& rule,
                    const Real& t)
        : FirstOrderReactionEvent(world, rule, t),
          product_(rule_.products().at(0))
    {
        // assert(rule_.products().size() == 1);
    }

    virtual boost::optional<ReactionInfo> react(const ReactionInfo::identified_voxel& voxel);

protected:

    const Species product_;

};

struct DesorptionEvent : FirstOrderReactionEvent
{
    DesorptionEvent(boost::shared_ptr<SpatiocyteWorld> world,
                    const ReactionRule& rule,
                    const Species& product,
                    const std::string& serial_location_product,
                    const Real& t)
        : FirstOrderReactionEvent(world, rule, t),
          product_(product),
          serial_location_product_(serial_location_product)
    {
        // assert(rule_.products().size() == 1);
    }

    virtual boost::optional<ReactionInfo> react(const ReactionInfo::identified_voxel& voxel);

protected:

    const Species product_;
    const std::string serial_location_product_;
};


struct EliminationEvent : FirstOrderReactionEvent
{
    EliminationEvent(boost::shared_ptr<SpatiocyteWorld> world,
                     const ReactionRule& rule,
                     const Species& product1,
                     const Species& product2,
                     const std::string& serial_location_product,
                     const Real& t)
        : FirstOrderReactionEvent(world, rule, t),
          product1_(product1),
          product2_(product2),
          serial_location_product_(serial_location_product)
    {
        // assert(rule_.products().size() == 2)
    }

    virtual boost::optional<ReactionInfo> react(const ReactionInfo::identified_voxel& voxel);

protected:

    const Species product1_;
    const Species product2_;
    const std::string serial_location_product_;

};

struct GenerationBesideEvent : FirstOrderReactionEvent
{
    GenerationBesideEvent(boost::shared_ptr<SpatiocyteWorld> world,
                          const ReactionRule& rule,
                          const Species& product,
                          const std::string& serial_location_product,
                          const Real& t)
        : FirstOrderReactionEvent(world, rule, t),
          product_(product),
          serial_location_product_(serial_location_product)
    {
        // assert(rule_.products().size() == 2)
    }

    virtual boost::optional<ReactionInfo> react(const ReactionInfo::identified_voxel& voxel);

protected:

    const Species product_;
    const std::string serial_location_product_;

};

struct GenerationTwoEvent : FirstOrderReactionEvent
{
    GenerationTwoEvent(boost::shared_ptr<SpatiocyteWorld> world,
                       const ReactionRule& rule,
                       const Species& product1,
                       const Species& product2,
                       const std::string& serial_location_product,
                       const Real& t)
        : FirstOrderReactionEvent(world, rule, t),
          product1_(product1),
          product2_(product2),
          serial_location_product_(serial_location_product)
    {
        // assert(rule_.products().size() == 2)
    }

    virtual boost::optional<ReactionInfo> react(const ReactionInfo::identified_voxel& voxel);

protected:

    const Species product1_;
    const Species product2_;
    const std::string serial_location_product_;

};

FirstOrderReactionEvent*
generate_first_order_reaction_event(boost::shared_ptr<SpatiocyteWorld> world,
                                    const ReactionRule& rule,
                                    const Real& t);

} // spatiocyte

} // ecell4

#endif /* ECELL4_SPATIOCYTE_FIRST_ORDER_REACTION_EVENT_HPP */
