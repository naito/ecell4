#ifndef ECELL4_SPATIOCYTE_REACTIONS_HPP
#define ECELL4_SPATIOCYTE_REACTIONS_HPP

#include <boost/shared_ptr.hpp>
#include <ecell4/core/VoxelPool.hpp>
#include <ecell4/core/ReactionRule.hpp>

namespace ecell4
{

namespace spatiocyte
{

class ReactionInfo
{
public:

    typedef std::pair<ParticleID, Voxel> identified_voxel;
    typedef std::vector<identified_voxel> container_type;

public:

    ReactionInfo() : t_(0), reactants_(), products_() {}

    ReactionInfo(const Real t) : t_(t), reactants_(), products_() {}

    ReactionInfo(
        const Real t,
        const container_type& reactants,
        const container_type& products)
        : t_(t), reactants_(reactants), products_(products) {}

    ReactionInfo(const ReactionInfo& another)
        : t_(another.t()), reactants_(another.reactants()), products_(another.products()) {}

    Real t() const
    {
        return t_;
    }

    bool has_occurred() const
    {
        return reactants_.size() > 0 || products_.size() > 0;
    }

    const container_type& reactants() const
    {
        return reactants_;
    }

    void add_reactant(const identified_voxel& pid_pair)
    {
        reactants_.push_back(pid_pair);
    }

    const container_type& products() const
    {
        return products_;
    }

    void add_product(const identified_voxel& pid_pair)
    {
        products_.push_back(pid_pair);
    }

protected:

    Real t_;
    container_type reactants_, products_;
};

} // spatiocyte

} // ecell4

#endif /* ECELL4_SPATIOCYTE_REACTIONS_HPP */
