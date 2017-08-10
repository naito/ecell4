#ifndef __ECELL4_SPATIOCYTE_UTILS_HPP
#define __ECELL4_SPATIOCYTE_UTILS_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>

namespace ecell4
{

namespace spatiocyte
{

template <typename T>
inline const T&
pick(const std::vector<T>& container, const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    assert(!container.empty());
    if (container.size() == 1)
        return container[0];
    return container[rng->uniform_int(0, container.size()-1)];
}

} // spatiocyte

} // ecell4

#endif /* __ECELL4_SPATIOCYTE_UTILS_HPP */
