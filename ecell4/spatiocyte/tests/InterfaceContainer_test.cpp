#define BOOST_TEST_MODULE "InterfaceContainer_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include "../InterfaceContainer.hpp"

using namespace ecell4;
using namespace ecell4::spatiocyte;

BOOST_AUTO_TEST_CASE(Constructor)
{
    InterfaceContainer container;
}

BOOST_AUTO_TEST_CASE(AddInterface)
{
    InterfaceContainer container;

    container.add_interface(0, 3);
    InterfaceContainer::iterator itr(container.find(0));

    BOOST_CHECK(!container.is_end(itr));
    BOOST_CHECK_EQUAL(1, container.num_adjoinings(itr));
    BOOST_CHECK_EQUAL(3, container.get_adjoining(itr, 0));

    container.add_interface(0, 5);
    BOOST_CHECK_EQUAL(2, container.num_adjoinings(itr));
    BOOST_CHECK_EQUAL(3, container.get_adjoining(itr, 0));
    BOOST_CHECK_EQUAL(5, container.get_adjoining(itr, 1));
}

BOOST_AUTO_TEST_CASE(Iterator)
{
    InterfaceContainer container;

    container.add_interface(0, 1);

    BOOST_CHECK(!container.is_end(container.find(0)));
    BOOST_CHECK(container.is_end(container.find(1)));
}
