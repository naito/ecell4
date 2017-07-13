#ifndef __ECELL4_SPATIOCYTE_INTERFACE_CONTAINER
#define __ECELL4_SPATIOCYTE_INTERFACE_CONTAINER

#include <vector>
#include <ecell4/core/get_mapper_mf.hpp>

namespace ecell4
{

namespace spatiocyte
{

template<typename coordinate_type>
class OneToManyMap
{

protected:
    typedef typename utils::get_mapper_mf<coordinate_type, std::vector<coordinate_type> >::type
            container_type;

    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

public:
    OneToManyMap() {}

    void add(coordinate_type interface, coordinate_type target)
    {
        iterator itr(container_.find(interface));

        if (itr != container_.end())
        {
            (*itr).second.push_back(target);
        }
        else
        {
            container_.insert(std::make_pair(interface, std::vector<coordinate_type>(1, target)));
        }
    }

    std::vector<coordinate_type> get(const coordinate_type& coordinate) const
    {
        const_iterator itr(container_.find(coordinate));

        if (itr == container_.end())
            return std::vector<coordinate_type>();
        else
            return (*itr).second;
    }

protected:
    container_type container_;

}; // class OneToManyMap

} // namespace spatiocyte

} // namespace ecell4

#endif /* __ECELL4_SPATIOCYTE_INTERFACE_CONTAINER */
