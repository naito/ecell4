#ifndef ECELL4_SPATIOCYTE_UTILS_HPP
#define ECELL4_SPATIOCYTE_UTILS_HPP

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

const Real calculate_dimensional_factor(
    const Shape::dimension_kind dimensionA,
    const Real D_A,
    bool is_structureA,
    const Shape::dimension_kind dimensionB,
    const Real D_B,
    bool is_structureB,
    boost::shared_ptr<const SpatiocyteWorld> world);

const Real calculate_alpha(
    const ReactionRule& rr,
    boost::shared_ptr<const SpatiocyteWorld> world);

} // spatiocyte

} // ecell4

#endif /* ECELL4_SPATIOCYTE_UTILS_HPP */
