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
