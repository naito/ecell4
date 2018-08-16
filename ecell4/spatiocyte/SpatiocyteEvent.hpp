#ifndef ECELL4_SPATIOCYTE_EVENT_HPP
#define ECELL4_SPATIOCYTE_EVENT_HPP

#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/Model.hpp>
#include "SpatiocyteReactions.hpp"
#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

struct SpatiocyteEvent : public Event
{
public:
    typedef std::pair<ReactionRule, ReactionInfo> reaction_type;

    SpatiocyteEvent(Real const& time) : Event(time) {}
    virtual ~SpatiocyteEvent() {}

    const std::vector<reaction_type>& reactions() const
    {
        return reactions_;
    }

    virtual void fire() {
        reactions_.clear();
        fire_();
    }

protected:
    virtual void fire_() = 0;

    void push_reaction(const ReactionRule& rule, const ReactionInfo& info)
    {
        reactions_.push_back(std::make_pair(rule, info));
    }

    std::vector<reaction_type> reactions_;

};

struct StepEvent : SpatiocyteEvent
{
    StepEvent(boost::shared_ptr<Model> model,
              boost::shared_ptr<SpatiocyteWorld> world,
              SpatiocyteWorld::space_type space,
              const Species& species,
              const Real& t,
              const Real alpha=1.0);
    virtual ~StepEvent() {}

    Species const& species() const
    {
        return mpool_->species();
    }

    Real const& alpha() const
    {
        return alpha_;
    }

    void fire_()
    {
        walk(alpha_);
        time_ += dt_;
    }

    virtual void walk(const Real& alpha) = 0;

protected:

    void attempt_reaction_(
        const ParticleID& pid,
        const Voxel& src,
        const Voxel& dst,
        const Real& alpha);

protected:

    boost::shared_ptr<Model> model_;
    boost::shared_ptr<SpatiocyteWorld> world_;
    SpatiocyteWorld::space_type space_;
    MoleculePool* mpool_;

    const Real alpha_;
};

struct StepEvent3D : StepEvent
{
    StepEvent3D(boost::shared_ptr<Model> model,
                boost::shared_ptr<SpatiocyteWorld> world,
                SpatiocyteWorld::space_type space,
                const Species& species,
                const Real& t,
                const Real alpha=1.0);

    void walk(const Real& alpha);
};

struct StepEvent2D : StepEvent
{
    StepEvent2D(boost::shared_ptr<Model> model,
                boost::shared_ptr<SpatiocyteWorld> world,
                SpatiocyteWorld::space_type space,
                const Species& species,
                const Real& t,
                const Real alpha=1.0);

    void walk(const Real& alpha);
};

struct ZerothOrderReactionEvent : SpatiocyteEvent
{
    ZerothOrderReactionEvent(
        boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t);

    virtual ~ZerothOrderReactionEvent() {}
    virtual void fire_();

    Real draw_dt();
    virtual void interrupt(Real const& t)
    {
        time_ = t + draw_dt();
    }

protected:

    boost::shared_ptr<SpatiocyteWorld> world_;
    ReactionRule rule_;
};

struct FirstOrderReactionEvent : SpatiocyteEvent
{
    FirstOrderReactionEvent(
        boost::shared_ptr<SpatiocyteWorld> world, const ReactionRule& rule, const Real& t);

    virtual ~FirstOrderReactionEvent() {}
    virtual void fire_();

    Real draw_dt();
    virtual void interrupt(Real const& t)
    {
        time_ = t + draw_dt();
    }

protected:

    ReactionInfo::Item choice()
    {
        const Species& species(rule_.reactants().at(0));
        const MoleculePool* mt(world_->find_molecule_pool(species));

        const Integer i(rng_.lock()->uniform_int(0, mt->size() - 1));
        const VoxelPool::coordinate_id_pair_type& info(mt->at(i));

        // TODO: Calling coordinate2voxel() is invalid
        return ReactionInfo::Item(info.pid, species, world_->coordinate2voxel(info.coordinate));
    }

    boost::shared_ptr<SpatiocyteWorld> world_;
    boost::weak_ptr<RandomNumberGenerator> rng_;
    ReactionRule rule_;
};

} // spatiocyte

} // ecell4

#endif /* ECELL4_SPATIOCYTE_EVENT_HPP */
