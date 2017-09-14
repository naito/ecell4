#ifndef ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/VoxelPool.hpp>
#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/get_mapper_mf.hpp>

#include "SpatiocyteWorld.hpp"
#include "SpatiocyteEvent.hpp"

namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteSimulator
    : public SimulatorBase<Model, SpatiocyteWorld>
{
public:

    typedef SimulatorBase<Model, SpatiocyteWorld> base_type;
    typedef SpatiocyteEvent::reaction_type reaction_type;
    typedef EventSchedulerBase<SpatiocyteEvent> scheduler_type;
    typedef utils::get_mapper_mf<Species, Real>::type alpha_map_type;

public:

    SpatiocyteSimulator(boost::shared_ptr<Model> model,
                        boost::shared_ptr<SpatiocyteWorld> world)
        : base_type(model, world)
    {
        initialize();
    }

    SpatiocyteSimulator(boost::shared_ptr<SpatiocyteWorld> world)
        : base_type(world)
    {
        initialize();
    }

    virtual Real dt() const
    {
        return dt_;
    }

    void initialize();
    void finalize();
    void step();
    bool step(const Real& upto);

    virtual bool check_reaction() const
    {
        return last_reactions().size() > 0;
    }

    const std::vector<SpatiocyteEvent::reaction_type>& last_reactions() const
    {
        return last_reactions_;
    }

protected:

    inline
    boost::shared_ptr<SpatiocyteEvent>
    create_step_event(const Species& species, const Real& t, const Real& alpha)
    {
        const MoleculePool *mpool(world_->find_molecule_pool(species).first);
        if (mpool->get_dimension() == Shape::THREE)
        {
            return boost::shared_ptr<SpatiocyteEvent>(
                    new StepEvent3D(model_, world_, species, t, alpha));
        }
        else
        {
            return boost::shared_ptr<SpatiocyteEvent>(
                    new StepEvent2D(model_, world_, species, t, alpha));
        }
    }

    inline
    boost::shared_ptr<SpatiocyteEvent>
    create_zeroth_order_reaction_event(const ReactionRule& reaction_rule, const Real& t)
    {
        return boost::shared_ptr<SpatiocyteEvent>(
                new ZerothOrderReactionEvent(world_, reaction_rule, t));
    }

    inline
    boost::shared_ptr<SpatiocyteEvent>
    create_first_order_reaction_event(const ReactionRule& reaction_rule, const Real& t)
    {
        return boost::shared_ptr<SpatiocyteEvent>(
                new FirstOrderReactionEvent(world_, reaction_rule, t));
    }

    void step_();
    void register_events(const Species& species);
    void update_alpha_map();

protected:

    scheduler_type scheduler_;
    alpha_map_type alpha_map_;

    std::vector<reaction_type> last_reactions_;

    Real dt_;
};

} // spatiocyte

} // ecell4

#endif /* ECELL4_LATTICE_LATTICE_SIMULATOR_HPP */
