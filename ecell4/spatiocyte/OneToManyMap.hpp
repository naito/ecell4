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

public:
    typedef typename container_type::iterator iterator;

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

    iterator find(const coordinate_type& coordinate)
    {
        return container_.find(coordinate);
    }

    bool is_end(iterator itr)
    {
        return itr == container_.end();
    }

    std::size_t num_adjoinings(iterator itr)
    {
        return (*itr).second.size();
    }

    coordinate_type get_adjoining(iterator itr, std::size_t idx)
    {
        return (*itr).second.at(idx);
    }

protected:
    container_type container_;

}; // class OneToManyMap

} // namespace spatiocyte

} // namespace ecell4

#endif /* __ECELL4_SPATIOCYTE_INTERFACE_CONTAINER */
