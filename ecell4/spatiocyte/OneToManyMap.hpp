#ifndef __ECELL4_SPATIOCYTE_INTERFACE_CONTAINER
#define __ECELL4_SPATIOCYTE_INTERFACE_CONTAINER

#include <vector>
#include <ecell4/core/get_mapper_mf.hpp>

namespace ecell4
{

namespace spatiocyte
{

template<typename T>
class OneToManyMap
{

protected:
    typedef typename utils::get_mapper_mf<T, std::vector<T> >::type
            container_type;

    typedef typename container_type::iterator iterator;

public:
    typedef typename container_type::const_iterator const_iterator;

    OneToManyMap() {}

    void add(T key, T value)
    {
        iterator itr(container_.find(key));

        if (itr != container_.end())
            (*itr).second.push_back(value);
        else
            container_.insert(std::make_pair(key, std::vector<T>(1, value)));
    }

    void add(T key, std::vector<T> values)
    {
        iterator itr(container_.find(key));

        if (itr != container_.end())
            std::copy(values.begin(), values.end(), back_inserter((*itr).second));
        else
            container_.insert(std::make_pair(key, values));
    }

    std::vector<T> get(const T& key) const
    {
        const_iterator itr(container_.find(key));

        if (itr == container_.end())
            return std::vector<T>();
        else
            return (*itr).second;
    }

    const_iterator begin() const { return container_.begin(); }
    const_iterator end() const { return container_.end(); }
    const_iterator find(const T& key) const { return container_.find(key); }

protected:
    container_type container_;

}; // class OneToManyMap

} // namespace spatiocyte

} // namespace ecell4

#endif /* __ECELL4_SPATIOCYTE_INTERFACE_CONTAINER */
