#ifndef ECELL4_VACANT_TYPE_HPP
#define ECELL4_VACANT_TYPE_HPP
#include "StructureType.hpp"

namespace ecell4
{

class VacantType
    : public StructureType
{
private:

    typedef StructureType base_type;

public:

    VacantType(const Shape::dimension_kind& dimension=Shape::THREE)
        : base_type(Species("", "0", "0"), NULL, 0, dimension)
    {
        ; // do nothing
    }

    ~VacantType()
    {
        ; // do nothing
    }

    voxel_type_type const voxel_type() const
    {
        return VACANT;
    }
};

} // ecell4

#endif /* ECELL4_VACANT_TYPE_HPP */
