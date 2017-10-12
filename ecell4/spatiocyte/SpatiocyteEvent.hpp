#ifndef ECELL4_SPATIOCYTE_EVENT_HPP
#define ECELL4_SPATIOCYTE_EVENT_HPP

#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/Model.hpp>
#include "ReactionInfo.hpp"
#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

struct SpatiocyteEvent : public Event
{
public:
    typedef std::pair<ReactionRule, ReactionInfo> reaction_type;

    SpatiocyteEvent(boost::shared_ptr<SpatiocyteWorld> world, Real const& time)
        : Event(time), world_(world) {}
    virtual ~SpatiocyteEvent() {}

    const std::vector<reaction_type>& reactions() const
    {
        return reactions_;
    }

    const std::vector<Species>& new_species() const
    {
        return new_species_;
    }

    virtual void fire() {
        reactions_.clear();
        new_species_.clear();
        fire_();
        filter_generated_species();
    }

    virtual void finalize(const Real& t) {}

protected:
    virtual void fire_() = 0;

    void push_reaction(const reaction_type& reaction)
    {
        reactions_.push_back(reaction);
    }

    void push_product(const Species& product)
    {
        if (!world_->has_species(product))
            new_species_.push_back(product);
    }

    boost::shared_ptr<SpatiocyteWorld> world_;

private:

    void filter_generated_species()
    {
        std::sort(new_species_.begin(), new_species_.end());
        std::vector<Species>::iterator end(std::unique(new_species_.begin(), new_species_.end()));
        new_species_.resize(std::distance(new_species_.begin(), end));

        std::vector<Species> filtered(new_species_.size());
        for (std::vector<Species>::const_iterator itr(new_species_.begin());
             itr != new_species_.end(); ++itr)
        {
            if (world_->has_species(*itr))
                filtered.push_back(*itr);
        }
        new_species_ = filtered;
    }

    std::vector<reaction_type> reactions_;
    std::vector<Species> new_species_;

};

struct StepEvent : SpatiocyteEvent
{
    StepEvent(boost::shared_ptr<Model> model,
              boost::shared_ptr<SpatiocyteWorld> world,
              const Species& species, const Real& t, const Real& alpha)
        : SpatiocyteEvent(world, t), model_(model), species_(species), alpha_(alpha) {}
    virtual ~StepEvent() {}

    virtual void fire_()
    {
        walk(alpha_);
        time_ += dt_;
    }

    virtual void finalize(const Real& t)
    {
        const Real queued_time(time() - dt());
        walk(alpha() * (t - queued_time) / dt());
    }

    Species const& species() const
    {
        return species_;
    }

    Real const& alpha() const
    {
        return alpha_;
    }

protected:

    virtual void walk(const Real& alpha) = 0;

    void attempt_reaction_(const SpatiocyteWorld::coordinate_id_pair_type& info,
                           const SpatiocyteWorld::coordinate_type to_coord,
                           const Real& alpha);

    boost::shared_ptr<Model> model_;
    Species species_;
    const Real alpha_;
};

struct StepEvent3D : StepEvent
{
protected:
    typedef StepEvent base;

public:
    StepEvent3D(boost::shared_ptr<Model> model,
                boost::shared_ptr<SpatiocyteWorld> world,
                const Species& species, const Real& t, const Real& alpha);

    virtual void walk(const Real& alpha);
};

struct StepEvent2D: StepEvent
{
protected:
    typedef StepEvent base;

public:
    StepEvent2D(boost::shared_ptr<Model> model,
                boost::shared_ptr<SpatiocyteWorld> world,
                const Species& species, const Real& t, const Real& alpha);

    virtual void walk(const Real& alpha);
};

struct ZerothOrderReactionEvent : SpatiocyteEvent
{
    ZerothOrderReactionEvent(
        boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t);

    virtual ~ZerothOrderReactionEvent() {}
    virtual void fire_();

    Real draw_dt();

protected:

    const ReactionRule rule_;
};


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

    virtual ReactionInfo react(const ReactionInfo::identified_voxel& voxel) = 0;

    virtual void fire_()
    {
        if (world_->num_voxels_exact(reactant_) != 0)
        {
            const ReactionInfo rinfo(react(world_->choice(reactant_)));
            if (rinfo.has_occurred())
            {
                push_reaction(std::make_pair(rule_, rinfo));
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

    virtual ~VanishmentEvent() {}
    virtual ReactionInfo react(const ReactionInfo::identified_voxel& voxel);
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

    virtual ~RearrangementEvent() {}
    virtual ReactionInfo react(const ReactionInfo::identified_voxel& voxel);

protected:

    const Species product_;

};

struct EliminationEvent : FirstOrderReactionEvent
{
    EliminationEvent(boost::shared_ptr<SpatiocyteWorld> world,
                     const ReactionRule& rule,
                     const Real& t)
        : FirstOrderReactionEvent(world, rule, t),
          product1_(rule_.products().at(0)),
          product2_(rule_.products().at(1))
    {
        // assert(rule_.products().size() == 2)
    }

    virtual ~EliminationEvent() {}
    virtual ReactionInfo react(const ReactionInfo::identified_voxel& voxel);

protected:

    const Species& product1_;
    const Species& product2_;

};

inline const std::string
get_serial(boost::shared_ptr<SpatiocyteWorld> world, const Integer coord)
{
    return world->get_voxel_pool_at(coord)->get_serial();
}

inline const std::string
get_location(boost::shared_ptr<SpatiocyteWorld> world, const Integer coord)
{
    return world->get_voxel_pool_at(coord)->get_location_serial();
}

} // spatiocyte

} // ecell4

#endif /* ECELL4_SPATIOCYTE_EVENT_HPP */
