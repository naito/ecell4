#define BOOST_TEST_MODULE "OneToManyMap_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include "../OneToManyMap.hpp"

using namespace ecell4;
using namespace ecell4::spatiocyte;

BOOST_AUTO_TEST_CASE(Constructor)
{
    OneToManyMap<int> container;
}

BOOST_AUTO_TEST_CASE(AddInterface)
{
    OneToManyMap<int> container;

    container.add(0, 3);
    std::vector<int> values(container.get(0));

    BOOST_CHECK(!values.empty());
    BOOST_CHECK_EQUAL(1, values.size());
    BOOST_CHECK_EQUAL(3, values.at(0));

    container.add(0, 5);
    BOOST_CHECK_EQUAL(1, values.size());

    values = container.get(0);
    BOOST_CHECK_EQUAL(2, values.size());
    BOOST_CHECK_EQUAL(3, values.at(0));
    BOOST_CHECK_EQUAL(5, values.at(1));
}

BOOST_AUTO_TEST_CASE(GET)
{
    OneToManyMap<int> container;

    container.add(0, 1);

    BOOST_CHECK(!container.get(0).empty());
    BOOST_CHECK(container.get(1).empty());
}
