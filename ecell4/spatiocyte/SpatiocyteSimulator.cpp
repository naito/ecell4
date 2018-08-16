#include "SpatiocyteSimulator.hpp"
#include "utils.hpp"

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
    const std::vector<std::pair<SpatiocyteWorld::space_type, Species> >
        space_species(world_->list_space_species());
    for (std::vector<std::pair<SpatiocyteWorld::space_type, Species> >::const_iterator
            itr(space_species.begin());
         itr != space_species.end(); ++itr)
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
        const boost::shared_ptr<SpatiocyteEvent>
            zeroth_order_reaction_event(
                create_zeroth_order_reaction_event(rr, world_->t()));
        scheduler_.add(zeroth_order_reaction_event);
    }

    dt_ = scheduler_.next_time() - t();
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

        const Real alpha(calculate_alpha(*itr, world_));
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

void SpatiocyteSimulator::register_events(
        const std::pair<SpatiocyteWorld::space_type, Species>& space_species_pair)
{
    const SpatiocyteWorld::space_type space(space_species_pair.first);
    const Species& sp(space_species_pair.second);

    const VoxelPool* pool(space->find_voxel_pool(sp));
    if (!pool->is_structure())
    {
        alpha_map_type::const_iterator itr(alpha_map_.find(sp));
        const Real alpha(itr != alpha_map_.end() ? itr->second : 1.0);
        scheduler_.add(create_step_event(space, sp, world_->t(), alpha));
    }

    const std::vector<ReactionRule>& reaction_rules(model_->query_reaction_rules(sp));
    for (std::vector<ReactionRule>::const_iterator itr(reaction_rules.begin());
        itr != reaction_rules.end(); ++itr)
    {
        scheduler_.add(create_first_order_reaction_event(*itr, world_->t()));
    }
}

boost::shared_ptr<SpatiocyteEvent> SpatiocyteSimulator::create_step_event(
        const SpatiocyteWorld::space_type space,
        const Species& species, const Real& t, const Real& alpha)
{
    const MoleculePool* mpool(world_->find_molecule_pool(species));

    if (mpool->get_dimension() == Shape::THREE)
    {
        return boost::shared_ptr<SpatiocyteEvent>(
                new StepEvent3D(model_, world_, space, species, t, alpha));
    }
    else if (mpool->get_dimension() == Shape::TWO)
    {
        return boost::shared_ptr<SpatiocyteEvent>(
                new StepEvent2D(model_, world_, space, species, t, alpha));
    }
    else
    {
        throw NotSupported("The dimension of a structure must be two or three.");
    }
}

boost::shared_ptr<SpatiocyteEvent>
SpatiocyteSimulator::create_zeroth_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<SpatiocyteEvent> event(
            new ZerothOrderReactionEvent(world_, reaction_rule, t));
    return event;
}

boost::shared_ptr<SpatiocyteEvent>
SpatiocyteSimulator::create_first_order_reaction_event(
    const ReactionRule& reaction_rule, const Real& t)
{
    boost::shared_ptr<SpatiocyteEvent> event(new FirstOrderReactionEvent(
                world_, reaction_rule, t));
    return event;
}

void SpatiocyteSimulator::finalize()
{
    scheduler_type::events_range events(scheduler_.events());
    for (scheduler_type::events_range::iterator itr(events.begin());
            itr != events.end(); ++itr)
    {
        const Real queued_time((*itr).second->time() - (*itr).second->dt());
        StepEvent* step_event(dynamic_cast<StepEvent*>((*itr).second.get()));
        if (step_event != NULL && queued_time < t())
        {
            const Real factor((t() - queued_time) / (*itr).second->dt());
            // assert(factor <= 1);
            step_event->walk(step_event->alpha() * factor);
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
    scheduler_type::value_type top(scheduler_.pop());
    const Real time(top.second->time());
    world_->set_t(time);
    world_->clear_new_species();
    top.second->fire(); // top.second->time_ is updated in fire()
    set_last_event_(boost::const_pointer_cast<const SpatiocyteEvent>(top.second));

    last_reactions_ = last_event_->reactions();

    scheduler_type::events_range events(scheduler_.events());
    for (scheduler_type::events_range::iterator itr(events.begin());
        itr != events.end(); ++itr)
    {
        (*itr).second->interrupt(time);
        scheduler_.update(*itr);
    }
    scheduler_.add(top.second);

    // update_alpha_map(); // may be performance cost
    const std::vector<std::pair<SpatiocyteWorld::space_type, Species> >&
        new_species(world_->get_new_species());
    for (std::vector<std::pair<SpatiocyteWorld::space_type, Species> >::const_iterator itr(new_species.begin());
         itr != new_species.end(); ++itr)
    {
        register_events(*itr);
    }

    num_steps_++;
}

} // spatiocyte

} // ecell4
