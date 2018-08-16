#include "utils.hpp"


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
    boost::shared_ptr<const SpatiocyteWorld> world)
{
    const Real voxel_radius(world->voxel_radius());
    const Real sqrt2(sqrt(2.0));

    if (dimensionA == Shape::THREE && dimensionB == Shape::THREE)
    {
        return 1. / (6 * sqrt2 * (D_A + D_B) * voxel_radius);
    }

    if (dimensionA == Shape::TWO && dimensionB == Shape::TWO)
    {
        const Real sqrt3(sqrt(3.0));
        const Real sqrt6(sqrt(6.0));
        const Real gamma(pow(2*sqrt2 + 4*sqrt3 + 3*sqrt6 + sqrt(22.0), 2) /
                        (72*(6*sqrt2 + 4*sqrt3 + 3*sqrt6)));
        return gamma / (D_A + D_B);
    }

    if (dimensionA == Shape::THREE && dimensionB == Shape::TWO)
    {
        Real factor = sqrt2 / (3 * D_A * voxel_radius);
        if (is_structureB)
        {
            factor *= world->unit_area();
        }
        return factor;
    }

    if (dimensionA == Shape::TWO && dimensionB == Shape::THREE)
    {
        Real factor = sqrt2 / (3 * D_B * voxel_radius);
        if (is_structureA)
        {
            factor *= world->unit_area();
        }
        return factor;
    }

    throw NotSupported("The dimension of a structure must be two or three.");
}

struct DimensionInfo
{
    DimensionInfo (Shape::dimension_kind dimension, Real D, bool is_structure)
        : dimension(dimension), D(D), is_structure(is_structure) {}

    Shape::dimension_kind dimension;
    Real D;
    bool is_structure;
};

inline DimensionInfo
get_dimension_info(
    boost::shared_ptr<const SpatiocyteWorld> world,
    const Species& species)
{
    try
    {
        const VoxelPool* pool = world->find_voxel_pool(species);
        return DimensionInfo(pool->get_dimension(), pool->D(), pool->is_structure());
    }
    catch (NotFound e) {}

    const MoleculeInfo info = world->get_molecule_info(species);

    if (info.loc != "")
    {
        try
        {
            const VoxelPool* location = world->find_voxel_pool(Species(info.loc));
            return DimensionInfo(location->get_dimension(), info.D, false);
        } catch (NotFound e) {}
    }

    return DimensionInfo(Shape::THREE, info.D, false);
}

const Real calculate_alpha(const ReactionRule& rr, boost::shared_ptr<const SpatiocyteWorld> world)
{
    const ReactionRule::reactant_container_type& reactants(rr.reactants());
    if (reactants.size() != 2)
        return 1.0;

    const DimensionInfo dinfoA(get_dimension_info(world, reactants.at(0)));
    const DimensionInfo dinfoB(get_dimension_info(world, reactants.at(1)));

    const Real factor(calculate_dimensional_factor(
                dinfoA.dimension, dinfoA.D, dinfoA.is_structure,
                dinfoB.dimension, dinfoB.D, dinfoB.is_structure,
                world));

    const Real alpha(1.0 / (factor * rr.k()));

    return alpha < 1.0 ? alpha : 1.0;
}

} // spatiocyte

} // ecell4
