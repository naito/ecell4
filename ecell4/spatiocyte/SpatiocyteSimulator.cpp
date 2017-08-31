#include "SpatiocyteSimulator.hpp"

#include <algorithm>
#include <iterator>
#include <ecell4/core/StructureType.hpp>

namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteSimulator::initialize()
{
    last_reactions_.clear();

    scheduler_.clear();
    update_alpha_map();

    const std::vector<Species> species(world_->list_species());
    for (std::vector<Species>::const_iterator itr(species.begin());
        itr != species.end(); ++itr)
    {
        register_events(*itr);
    }

    const std::vector<ReactionRule>& rules(model_->reaction_rules());
    for (std::vector<ReactionRule>::const_iterator i(rules.begin());
        i != rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        if (rr.reactants().size() != 0)
        {
            continue;
        }
        scheduler_.add(create_zeroth_order_reaction_event(rr, world_->t()));
    }

    dt_ = scheduler_.next_time() - t();
}

inline const Real
calculate_alpha(const ReactionRule::reactant_container_type& reactants,
                const Real& k,
                const boost::shared_ptr<SpatiocyteWorld>& world)
{
    try
    {
        const VoxelPool *vp0(world->find_voxel_pool(reactants.at(0)));
        const VoxelPool *vp1(world->find_voxel_pool(reactants.at(1)));

        const Real dfactor(calculate_dimensional_factor(vp0, vp1,
                           boost::const_pointer_cast<const SpatiocyteWorld>(world)));
        const Real inv_alpha = dfactor * k;

        return inv_alpha <= 1.0 ? 1.0 : 1.0 / inv_alpha;
    }
    catch(NotFound e)
    {
        return 1.0;
    }
}

void SpatiocyteSimulator::update_alpha_map()
{
    boost::shared_ptr<Model> model_(model());
    if (!model_ || !model_->is_static())
        return;

    const Model::reaction_rule_container_type reaction_rules(model_->reaction_rules());
    for (Model::reaction_rule_container_type::const_iterator itr(reaction_rules.begin());
            itr != reaction_rules.end(); ++itr)
    {
        const ReactionRule::reactant_container_type& reactants((*itr).reactants());
        if (reactants.size() != 2)
            continue;

        const Real alpha(calculate_alpha(reactants, (*itr).k(), world_));

        for (int i(0); i < 2; ++i) {
            const Species& sp(reactants.at(i));
            alpha_map_type::iterator map_itr(alpha_map_.find(sp));

            if (map_itr == alpha_map_.end())
                alpha_map_.insert(alpha_map_type::value_type(sp, alpha));
            else if ((*map_itr).second > alpha)
                (*map_itr).second = alpha;
        }
    }
}

void SpatiocyteSimulator::register_events(const Species& sp)
{
    if (world_->has_molecule_pool(sp))
    {
        //TODO: Call steps only if sp is assigned not to StructureType.
        alpha_map_type::const_iterator itr(alpha_map_.find(sp));
        const Real alpha(itr != alpha_map_.end() ? itr->second : 1.0);
        scheduler_.add(create_step_event(sp, world_->t(), alpha));
    }

    std::vector<ReactionRule> reaction_rules(model_->query_reaction_rules(sp));
    for (std::vector<ReactionRule>::const_iterator rr(reaction_rules.begin());
        rr != reaction_rules.end(); ++rr)
    {
        scheduler_.add(create_first_order_reaction_event(*rr, world_->t()));
    }
}

void SpatiocyteSimulator::finalize()
{
    scheduler_type::events_range events(scheduler_.events());
    for (scheduler_type::events_range::iterator event(events.begin());
            event != events.end(); ++event)
    {
        const Real queued_time((*event).second->time() - (*event).second->dt());
        StepEvent* step_event(dynamic_cast<StepEvent*>((*event).second.get()));
        if (step_event != NULL && queued_time < t())
        {
            const Real alpha((t() - queued_time) / (*event).second->dt());
            step_event->walk(step_event->alpha() * alpha);
        }
    }

    initialize();
}

void SpatiocyteSimulator::step()
{
    step_();
    dt_ = scheduler_.next_time() - t();
}

bool SpatiocyteSimulator::step(const Real& upto)
{
    if (upto < t())
    {
        return false;
    }

    if (scheduler_.size() > 0 && upto >= scheduler_.top().second->time())
    {
        step_();
        dt_ = scheduler_.next_time() - t();
        return true;
    }

    world_->set_t(upto);
    last_reactions_.clear();
    dt_ = scheduler_.next_time() - t();
    finalize();
    return false;
}

void SpatiocyteSimulator::step_()
{
    boost::shared_ptr<SpatiocyteEvent> next_event(scheduler_.pop().second);
    const Real time(next_event->time());
    world_->set_t(time);
    next_event->fire();
    // next_event->time_ is updated in fire()
    last_reactions_ = next_event->reactions();

    scheduler_type::events_range events(scheduler_.events());
    for (scheduler_type::events_range::iterator event(events.begin());
         event != events.end(); ++event)
    {
        (*event).second->interrupt(time);
        scheduler_.update(*event);
    }

    // update_alpha_map(); // may be performance cost
    const std::vector<Species>& new_species(next_event->new_species());
    for (std::vector<Species>::const_iterator species(new_species.begin());
        species != new_species.end(); ++species)
    {
        register_events(*species);
    }

    scheduler_.add(next_event);

    num_steps_++;
}

} // spatiocyte

} // ecell4
