#include "SpatiocyteEvent.hpp"
#include "utils.hpp"

namespace ecell4
{

namespace spatiocyte
{

StepEvent::StepEvent(boost::shared_ptr<Model> model, boost::shared_ptr<SpatiocyteWorld> world,
        const Species& species, const Real& t, const Real alpha)
    : SpatiocyteEvent(t),
      model_(model),
      world_(world),
      mpool_(world_->find_molecule_pool(species)),
      alpha_(alpha)
{
    time_ = t;
}

StepEvent3D::StepEvent3D(boost::shared_ptr<Model> model,
                         boost::shared_ptr<SpatiocyteWorld> world,
                         const Species& species,
                         const Real& t,
                         const Real alpha)
    : StepEvent(model, world, species, t, alpha)
{
    const SpatiocyteWorld::molecule_info_type minfo(world_->get_molecule_info(species));
    const Real D(minfo.D);
    const Real R(world_->voxel_radius());

    if (D <= 0)
        dt_ = inf;
    else
        dt_ = 2 * R * R / 3 / D * alpha_;

    time_ = t + dt_;
}

void StepEvent3D::walk(const Real& alpha)
{
    if (alpha < 0 || alpha > 1)
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MoleculePool::container_type voxels;
    copy(mpool_->begin(), mpool_->end(), back_inserter(voxels));

    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(voxels.begin());
         itr != voxels.end(); ++itr)
    {
        const SpatiocyteWorld::coordinate_id_pair_type& info(*itr);
        const Voxel voxel(world_->coordinate2voxel(info.coordinate));
        const Integer rnd(rng->uniform_int(0, voxel.num_neighbors()-1));

        // should skip if a voxel is not the target species.
        // when reaction has occured before, a voxel can be changed.
        if (voxel.get_voxel_pool() != mpool_)
        {
            continue;
        }

        const Voxel neighbor(voxel.get_neighbor(rnd));

        if (boost::optional<const std::vector<Voxel>&> targets = world_->find_interface(neighbor))
        {
            const Voxel new_neighbor(targets->at(rng->uniform_int(0, targets->size()-1)));
            attempt_reaction_(info.pid, voxel, new_neighbor, alpha);
        }
        else if (world_->can_move(voxel, neighbor))
        {
            if (rng->uniform(0,1) <= alpha)
            {
                world_->move(voxel, neighbor, /*candidate=*/idx);
            }
        }
        else
        {
            attempt_reaction_(info.pid, voxel, neighbor, alpha);
        }

        ++idx;
    }
}

StepEvent2D::StepEvent2D(boost::shared_ptr<Model> model,
                         boost::shared_ptr<SpatiocyteWorld> world,
                         const Species& species,
                         const Real& t,
                         const Real alpha)
    : StepEvent(model, world, species, t, alpha)
{
    const SpatiocyteWorld::molecule_info_type minfo(world_->get_molecule_info(species));
    const Real D(minfo.D);
    const Real R(world_->voxel_radius());

    if (D <= 0)
        dt_ = inf;
    else
        dt_ = R * R / D * alpha_;

    time_ = t + dt_;
}

void StepEvent2D::walk(const Real& alpha)
{
    if (alpha < 0 || alpha > 1)
    {
        return; // INVALID ALPHA VALUE
    }

    const boost::shared_ptr<RandomNumberGenerator>& rng(world_->rng());
    MoleculePool::container_type voxels;
    copy(mpool_->begin(), mpool_->end(), back_inserter(voxels));

    std::size_t idx(0);
    for (MoleculePool::container_type::iterator itr(voxels.begin());
         itr != voxels.end(); ++itr)
    {
        const SpatiocyteWorld::coordinate_id_pair_type& info(*itr);

        // TODO: Calling coordinate2voxel is invalid
        const Voxel voxel(world_->coordinate2voxel(info.coordinate));

        // should skip if a voxel is not the target species.
        // when reaction has occured before, a voxel can be changed.
        if (voxel.get_voxel_pool() != mpool_)
        {
            continue;
        }

        std::vector<std::size_t> nids;
        for (std::size_t i(0); i < voxel.num_neighbors(); ++i)
        {
            nids.push_back(i);
        }
        ecell4::shuffle(*(rng.get()), nids);

        for (std::vector<std::size_t>::const_iterator itr(nids.begin());
             itr != nids.end(); ++itr)
        {
            const Voxel neighbor(voxel.get_neighbor(*itr));
            boost::shared_ptr<const VoxelPool> target(neighbor.get_voxel_pool());

            if (target->get_dimension() > mpool_->get_dimension())
                continue;

            if (boost::optional<const std::vector<Voxel>&> targets = world_->find_interface(neighbor))
            {
                const Voxel new_neighbor(targets->at(rng->uniform_int(0, targets->size()-1)));
                attempt_reaction_(info.pid, voxel, new_neighbor, alpha);
            }
            else if (world_->can_move(voxel, neighbor))
            {
                if (rng->uniform(0,1) <= alpha)
                    world_->move(voxel, neighbor, /*candidate=*/idx);
            }
            else
            {
                attempt_reaction_(info.pid, voxel, neighbor, alpha);
            }
            break;
        }
        ++idx;
    }
}

void StepEvent::attempt_reaction_(
    const ParticleID& pid,
    const Voxel& src,
    const Voxel& dst,
    const Real& alpha)
{
    boost::shared_ptr<const VoxelPool> from_mt(src.get_voxel_pool());
    boost::shared_ptr<const VoxelPool> to_mt(dst.get_voxel_pool());

    if (to_mt->is_vacant())
    {
        return;
    }

    const Species& speciesA(from_mt->species());
    const Species& speciesB(to_mt->species());

    const std::vector<ReactionRule> rules(model_->query_reaction_rules(speciesA, speciesB));

    if (rules.empty())
    {
        return;
    }

    const Real factor(calculate_dimensional_factor(from_mt, to_mt, world_));
    const Real rnd(world_->rng()->uniform(0,1));
    Real accp(0.0);

    for (std::vector<ReactionRule>::const_iterator itr(rules.begin()); itr != rules.end(); ++itr)
    {
        const Real k((*itr).k());
        const Real P(k * factor * alpha);
        accp += P;
        if (accp > 1)
        {
            std::cerr << "The total acceptance probability [" << accp
                << "] exceeds 1 for '" << speciesA.serial()
                << "' and '" << speciesB.serial() << "'." << std::endl;
        }
        if (accp >= rnd)
        {
            ReactionInfo rinfo(apply_second_order_reaction(
                        world_, *itr,
                        ReactionInfo::Item(pid, from_mt->species(), src),
                        ReactionInfo::Item(to_mt->get_particle_id(dst.coordinate),
                                           to_mt->species(), dst)));
            if (rinfo.has_occurred())
            {
                reaction_type reaction(std::make_pair(*itr, rinfo));
                push_reaction(reaction);
            }
            return;
        }
    }
}

} // spatiocyte

} // ecell4
