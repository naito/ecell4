#ifndef ECELL4_SPATIOCYTE_VOXEL_HPP
#define ECELL4_SPATIOCYTE_VOXEL_HPP

#include <ecell4/core/hash.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/VoxelSpaceBase.hpp>
#include <boost/weak_ptr.hpp>

namespace ecell4
{

namespace spatiocyte
{

struct Voxel
{
    typedef Integer coordinate_type;
    typedef VoxelSpaceBase space_type;

    Voxel(boost::weak_ptr<space_type> space, coordinate_type coordinate)
        : space(space), coordinate(coordinate)
    {
    }

    boost::weak_ptr<space_type> space;
    coordinate_type coordinate;

public:

    const Real3 position() const
    { return space.lock()->coordinate2position(coordinate);
    }

    bool clear() const
    {
        return space.lock()->remove_voxel(coordinate);
    }

    boost::shared_ptr<VoxelPool> get_voxel_pool() const
    {
        return space.lock()->get_voxel_pool_at(coordinate);
    }

    Integer num_neighbors() const
    {
        return space.lock()->num_neighbors(coordinate);
    }

    Voxel get_neighbor(Integer nrand) const
    {
        return Voxel(space, space.lock()->get_neighbor(coordinate, nrand));
    }

    bool operator==(const Voxel& rhs) const
    {
        return this->space.lock() == rhs.space.lock() &&
               this->coordinate == rhs.coordinate;
    }
};

}

}

ECELL4_DEFINE_HASH_BEGIN()

template<>
struct hash<ecell4::spatiocyte::Voxel>
{
    typedef ecell4::spatiocyte::Voxel voxel_type;

    std::size_t operator()(const voxel_type& val) const
    {
        return hash<voxel_type::space_type*>()(val.space.lock().get()) ^
               hash<voxel_type::coordinate_type>()(val.coordinate);
    }
};

ECELL4_DEFINE_HASH_END()

#endif
