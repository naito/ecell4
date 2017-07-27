#ifndef __ECELL4_SPATIOCYTE_UTILS_HPP
#define __ECELL4_SPATIOCYTE_UTILS_HPP

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

inline const Real
calculate_dimensional_factor(const VoxelPool* vp_a,
                             const VoxelPool* vp_b,
                             const boost::shared_ptr<const SpatiocyteWorld>& world)
{
    const Real voxel_radius(world->voxel_radius());
    const Real D_a(vp_a->D()),
               D_b(vp_b->D());
    const Shape::dimension_kind dim_a(vp_a->get_dimension()),
                                dim_b(vp_b->get_dimension());

    if (dim_a == Shape::THREE && dim_b == Shape::THREE)
    {
        return 2. / (3. * world->unit_voxel_volume() * (D_a + D_b) * voxel_radius);
    }
    else if (dim_a == Shape::TWO && dim_b == Shape::TWO)
    {
        const Real gamma(pow(2 * sqrt(2.0) + 4 * sqrt(3.0) + 3 * sqrt(6.0) + sqrt(22.0), 2) /
                         (72 * (6 * sqrt(2.0) + 4 * sqrt(3.0) + 3 * sqrt(6.0))));
        return gamma / (D_a + D_b);
    }
    else if (dim_a == Shape::THREE && dim_b == Shape::TWO)
    {
        const Real factor(sqrt(2.0) / (3 * D_a * voxel_radius));
        if (vp_b->is_structure())
            return factor * world->unit_area();
        return factor;
    }
    else if (dim_a == Shape::TWO && dim_b == Shape::THREE)
    {
        const Real factor(sqrt(2.0) / (3 * D_b * voxel_radius));
        if (vp_a->is_structure())
            return factor * world->unit_area();
        return factor;
    }
    throw NotSupported("The dimension of a structure must be two or three.");
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

} // spatiocyte

} // ecell4

#endif /* __ECELL4_SPATIOCYTE_UTILS_HPP */
