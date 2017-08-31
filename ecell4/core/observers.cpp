#include "observers.hpp"


namespace ecell4
{

const Real Observer::next_time() const
{
    return inf;
}

void Observer::initialize(const boost::shared_ptr<Space>& space)
{
    ;
}

void Observer::finalize(const boost::shared_ptr<Space>& space)
{
    ;
}

void Observer::reset()
{
    num_steps_ = 0;
}

bool Observer::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    ++num_steps_;
    return true;
}

const Integer Observer::num_steps() const
{
    return num_steps_;
}

const Real FixedIntervalObserver::next_time() const
{
    return t0_ + dt_ * count_;
}

const Integer FixedIntervalObserver::count() const
{
    return count_;
}

void FixedIntervalObserver::initialize(const boost::shared_ptr<Space>& space)
{
    base_type::initialize(space);

    if (dt_ <= 0.0)
    {
        throw std::invalid_argument(
            "A step interval of FixedIntervalObserver must be positive.");
    }

    if (count_ == 0)
    {
        t0_ = space->t();
    }
    else
    {
        while (next_time() < space->t())
        {
            ++count_;
        }
    }
}

bool FixedIntervalObserver::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    ++count_;
    return base_type::fire(sim, space);
}

void FixedIntervalObserver::reset()
{
    base_type::reset();
    count_ = 0;
    t0_ = 0.0; //DUMMY
}

void NumberLogger::log(const boost::shared_ptr<Space>& space)
{
    data_container_type::value_type tmp;
    tmp.push_back(space->t());
    for (species_container_type::const_iterator i(targets.begin());
        i != targets.end(); ++i)
    {
        tmp.push_back(space->get_value(*i));
        // tmp.push_back(space->num_molecules(*i));
    }
    data.push_back(tmp);
}

void NumberLogger::save(const std::string& filename) const
{
    if (!is_directory(filename))
    {
        throw NotFound("The output path does not exists.");
    }

    std::ofstream ofs(filename.c_str(), std::ios::out);
    ofs << std::setprecision(17);

    for (species_container_type::const_iterator i(targets.begin());
         i != targets.end(); ++i)
    {
        ofs << ",\"" << (*i).serial() << "\"";
    }
    ofs << std::endl;

    for (data_container_type::const_iterator i(data.begin());
         i != data.end(); ++i)
    {
        std::vector<Real>::const_iterator j((*i).begin());
        ofs << (*j);
        ++j;

        for (; j != (*i).end(); ++j)
        {
            ofs << "," << (*j);
        }
        ofs << std::endl;
    }

    ofs.close();
}

void FixedIntervalNumberObserver::initialize(const boost::shared_ptr<Space>& space)
{
    base_type::initialize(space);
    logger_.initialize();
}

bool FixedIntervalNumberObserver::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    logger_.log(space);
    return base_type::fire(sim, space);
}

void FixedIntervalNumberObserver::reset()
{
    logger_.reset();
    base_type::reset();
}

NumberLogger::data_container_type FixedIntervalNumberObserver::data() const
{
    return logger_.data;
}

NumberLogger::species_container_type FixedIntervalNumberObserver::targets() const
{
    return logger_.targets;
}

void NumberObserver::initialize(const boost::shared_ptr<Space>& space)
{
    base_type::initialize(space);
    logger_.initialize();
    logger_.log(space);
}

void NumberObserver::finalize(const boost::shared_ptr<Space>& space)
{
    if (logger_.data.size() == 0 || logger_.data.back()[0] != space->t())
    {
        logger_.log(space);
    }
    base_type::finalize(space);
}

bool NumberObserver::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    if (sim->check_reaction())
    {
        logger_.log(space);
        return base_type::fire(sim, space);
    }
    return true;
}

void NumberObserver::reset()
{
    logger_.reset();
    base_type::reset();
}

NumberLogger::data_container_type NumberObserver::data() const
{
    return logger_.data;
}

NumberLogger::species_container_type NumberObserver::targets() const
{
    return logger_.targets;
}

const Real TimingObserver::next_time() const
{
    if (count_ >= static_cast<Integer>(t_.size()))
    {
        return inf;
    }
    return t_[count_];
}

void TimingObserver::initialize(const boost::shared_ptr<Space>& space)
{
    while (next_time() < space->t())
    {
        ++count_;
    }
}

bool TimingObserver::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    ++num_steps_;
    ++count_;
    return true;
}

void TimingObserver::reset()
{
    num_steps_ = 0;
    count_ = 0;
}

void TimingNumberObserver::initialize(const boost::shared_ptr<Space>& space)
{
    base_type::initialize(space);
    logger_.initialize();
}

bool TimingNumberObserver::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    logger_.log(space);
    return base_type::fire(sim, space);
}

void TimingNumberObserver::reset()
{
    logger_.reset();
    base_type::reset();
}

NumberLogger::data_container_type TimingNumberObserver::data() const
{
    return logger_.data;
}

NumberLogger::species_container_type TimingNumberObserver::targets() const
{
    return logger_.targets;
}

void FixedIntervalHDF5Observer::initialize(const boost::shared_ptr<Space>& space)
{
    base_type::initialize(space);
}

bool FixedIntervalHDF5Observer::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    if (!is_directory(filename()))
    {
        throw NotFound("The output path does not exists.");
    }

    space->save(filename());

    return base_type::fire(sim, space);
}

const std::string FixedIntervalHDF5Observer::filename(const Integer idx) const
{
    boost::format fmt(prefix_);

    if (fmt.expected_args() == 0)
    {
        return fmt.str();
    }
    else
    {
        return (fmt % idx).str();
    }
}

void FixedIntervalCSVObserver::initialize(const boost::shared_ptr<Space>& space)
{
    base_type::initialize(space);
    logger_.initialize();
}

bool FixedIntervalCSVObserver::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    log(space);
    return base_type::fire(sim, space);
}

void FixedIntervalCSVObserver::log(const boost::shared_ptr<Space>& space)
{
    if (!is_directory(filename()))
    {
        throw NotFound("The output path does not exists.");
    }

    std::ofstream ofs(filename().c_str(), std::ios::out);
    logger_.save(ofs, space);
    ofs.close();
}

const std::string FixedIntervalCSVObserver::filename() const
{
    boost::format fmt(prefix_);

    if (fmt.expected_args() == 0)
    {
        return fmt.str();
    }
    else
    {
        return (fmt % num_steps()).str();
    }
}

void FixedIntervalCSVObserver::reset()
{
    logger_.reset();
    base_type::reset();
}

void CSVObserver::initialize(const boost::shared_ptr<Space>& space)
{
    base_type::initialize(space);
    logger_.initialize();
    log(space);
}

bool CSVObserver::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    const bool retval = base_type::fire(sim, space); // Increment num_steps_ first.
    log(space);
    return retval;
}

void CSVObserver::log(const boost::shared_ptr<Space>& space)
{
    if (!is_directory(filename()))
    {
        throw NotFound("The output path does not exists.");
    }

    std::ofstream ofs(filename().c_str(), std::ios::out);
    logger_.save(ofs, space);
    ofs.close();
}

const std::string CSVObserver::filename() const
{
    boost::format fmt(prefix_);

    if (fmt.expected_args() == 0)
    {
        return fmt.str();
    }
    else
    {
        return (fmt % num_steps()).str();
    }
}

void CSVObserver::reset()
{
    logger_.reset();
    base_type::reset();
}

void TimeoutObserver::initialize(const boost::shared_ptr<Space>& space)
{
    base_type::initialize(space);
    duration_ = 0.0;
    time(&tstart_);
}

void TimeoutObserver::finalize(const boost::shared_ptr<Space>& space)
{
    base_type::finalize(space);
    acc_ += duration_;
}

bool TimeoutObserver::fire(const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    time_t tnow;
    time(&tnow);
    duration_ = difftime(tnow, tstart_);
    if (duration_ >= interval_)
    {
        return false;
    }
    return true;
}

void TimeoutObserver::reset()
{
    base_type::reset();
    duration_ = 0.0;
    acc_ = 0.0;
    time(&tstart_);
}

const Real FixedIntervalTrackingObserver::next_time() const
{
    return std::min(event_.next_time(), subevent_.next_time());
}

const Integer FixedIntervalTrackingObserver::num_steps() const
{
    return event_.num_steps + subevent_.num_steps;
}

const Integer FixedIntervalTrackingObserver::count() const
{
    return event_.count;
}

void FixedIntervalTrackingObserver::initialize(const boost::shared_ptr<Space>& space)
{
    event_.initialize(space->t());
    subevent_.initialize(space->t());

    if (pids_.size() == 0)
    {
        typedef std::vector<std::pair<ParticleID, Particle> > particle_id_pairs;
        for (std::vector<Species>::const_iterator i(species_.begin());
             i != species_.end(); ++i)
        {
            const Species& sp(*i);
            particle_id_pairs const particles(space->list_particles_exact(sp));
            pids_.reserve(pids_.size() + particles.size());
            for (particle_id_pairs::const_iterator j(particles.begin());
                j != particles.end(); ++j)
            {
                pids_.push_back((*j).first);
            }
        }

        prev_positions_.clear();
        prev_positions_.resize(pids_.size(), Real3(0, 0, 0));
        trajectories_.clear();
        trajectories_.resize(pids_.size(), std::vector<Real3>());
        strides_.clear();
        strides_.resize(pids_.size(), Real3(0, 0, 0));
        t_.clear();
    }
}

bool FixedIntervalTrackingObserver::fire(
    const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    if (subevent_.next_time() <= event_.next_time())
    {
        fire_subevent(sim, space);
    }
    else
    {
        fire_event(sim, space);
    }
    return true;
}

void FixedIntervalTrackingObserver::fire_subevent(
    const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    typedef std::vector<std::pair<ParticleID, Particle> > particle_id_pairs;

    const Real3& edge_lengths(space->edge_lengths());

    std::vector<Real3>::iterator j(prev_positions_.begin());
    std::vector<Real3>::iterator k(strides_.begin());
    for (std::vector<ParticleID>::iterator i(pids_.begin());
        i != pids_.end(); ++i, ++j, ++k)
    {
        if ((*i) == ParticleID() || space->has_particle(*i))
        {
            continue;
        }

        const Real3 pos((*j) - (*k));
        Real Lmin(threshold_);
        ParticleID newpid;

        for (std::vector<Species>::const_iterator l(species_.begin());
             l != species_.end(); ++l)
        {
            const Species& sp(*l);
            particle_id_pairs const particles(space->list_particles_exact(sp));
            for (particle_id_pairs::const_iterator m(particles.begin());
                m != particles.end(); ++m)
            {
                if (std::find(pids_.begin(), pids_.end(), (*m).first) == pids_.end())
                {
                    const Real L(distance(pos, (*m).second.position(), edge_lengths));
                    if (L < Lmin)
                    {
                        Lmin = L;
                        newpid = (*m).first;
                    }
                }
            }
        }

        (*i) = newpid;
    }

    if (resolve_boundary_)
    {
        const Real3 edge_lengths(space->actual_lengths());
        std::vector<Real3>::iterator j(prev_positions_.begin());
        std::vector<Real3>::iterator k(strides_.begin());
        for (std::vector<ParticleID>::const_iterator i(pids_.begin());
            i != pids_.end(); ++i, ++j, ++k)
        {
            if ((*i) != ParticleID() && space->has_particle(*i))
            {
                Real3& stride(*k);
                Real3 pos(stride + space->get_particle(*i).second.position());
                if (subevent_.num_steps > 0)
                {
                    const Real3& prev(*j);
                    for (unsigned int dim(0); dim != 3; ++dim)
                    {
                        const Real L(edge_lengths[dim]);
                        if (pos[dim] - prev[dim] >= L * 0.5)
                        {
                            stride[dim] -= L;
                            pos[dim] -= L;
                        }
                        else if (pos[dim] - prev[dim] <= L * -0.5)
                        {
                            stride[dim] += L;
                            pos[dim] += L;
                        }
                    }
                }
                (*j) = pos;
            }
        }
    }

    subevent_.fire();
}

void FixedIntervalTrackingObserver::fire_event(
    const Simulator* sim, const boost::shared_ptr<Space>& space)
{
    t_.push_back(space->t());

    const Real3 edge_lengths(space->actual_lengths());
    std::vector<Real3>::const_iterator j(prev_positions_.begin());
    std::vector<Real3>::const_iterator k(strides_.begin());
    std::vector<std::vector<Real3> >::iterator l(trajectories_.begin());
    for (std::vector<ParticleID>::const_iterator i(pids_.begin());
        i != pids_.end(); ++i)
    {
        if (space->has_particle(*i))
        {
            const Real3& stride(*k);
            Real3 pos(stride + space->get_particle(*i).second.position());

            if (resolve_boundary_ && subevent_.num_steps > 0)
            {
                const Real3& prev(*j);

                for (unsigned int dim(0); dim != 3; ++dim)
                {
                    const Real L(edge_lengths[dim]);
                    if (pos[dim] - prev[dim] >= L * 0.5)
                    {
                        pos[dim] -= L;
                    }
                    else if (pos[dim] - prev[dim] <= L * -0.5)
                    {
                        pos[dim] += L;
                    }
                }
            }

            (*l).push_back(pos);
        }
        ++j;
        ++k;
        ++l;
    }

    event_.fire();
}

void FixedIntervalTrackingObserver::reset()
{
    event_.reset();
    subevent_.reset();

    prev_positions_.clear();
    prev_positions_.resize(pids_.size(), Real3(0, 0, 0));
    trajectories_.clear();
    trajectories_.resize(pids_.size(), std::vector<Real3>());
    strides_.clear();
    strides_.resize(pids_.size(), Real3(0, 0, 0));
    t_.clear();
}

const std::vector<std::vector<Real3> >& FixedIntervalTrackingObserver::data() const
{
    return trajectories_;
}

const Integer FixedIntervalTrackingObserver::num_tracers() const
{
    return pids_.size();
}

const std::vector<Real>& FixedIntervalTrackingObserver::t() const
{
    return t_;
}

} // ecell4
