#define BOOST_TEST_MODULE "SpatiocyteWorld_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>
#include <fstream>

#include <ecell4/core/Sphere.hpp>
#include "../SpatiocyteWorld.hpp"

using namespace ecell4;
using namespace ecell4::spatiocyte;

typedef SpatiocyteWorld::coordinate_type coordinate_type;

const Real DEFAULT_VOXEL_RADIUS = 1e-8;

struct SpatiocyteWorldTestFixture
{
    const Real VOXEL_RADIUS;
    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    SpatiocyteWorld world;

    SpatiocyteWorldTestFixture()
        : VOXEL_RADIUS(1e-8),
          rng(new GSLRandomNumberGenerator()),
          world(Real3(1e-6, 1e-6, 1e-6), VOXEL_RADIUS, rng)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(SpatiocyteWorldTest, SpatiocyteWorldTestFixture)


BOOST_AUTO_TEST_CASE(Constructor) {}

BOOST_AUTO_TEST_CASE(DefaultValues)
{
    BOOST_CHECK_EQUAL(0, world.t());
    BOOST_CHECK_EQUAL(0, world.num_species());
    BOOST_CHECK_EQUAL(0, world.num_particles());
    BOOST_CHECK_EQUAL(0, world.num_voxels());
    BOOST_CHECK_EQUAL(0, world.list_particles().size());
}

BOOST_AUTO_TEST_CASE(Time)
{
    world.set_t(23.4);
    BOOST_CHECK_EQUAL(23.4, world.t());
}

BOOST_AUTO_TEST_CASE(HasSpecies)
{
    BOOST_CHECK(!world.has_species(Species("Species")));
}

BOOST_AUTO_TEST_CASE(UpdateParticle)
{
    SerialIDGenerator<ParticleID> sidgen;
    ParticleID pid(sidgen());
    Species species("A");

    BOOST_ASSERT(world.update_particle(pid,
                Particle(species,
                         /* position = */ Real3(2.0e-7, 1.0e-7, 0.0),
                         /*   radius = */ VOXEL_RADIUS,
                         /*        D = */ 1.0e-12)));

    BOOST_CHECK(world.has_species(species));
    BOOST_CHECK(world.has_particle(pid));
    BOOST_CHECK_EQUAL(1, world.list_particles().size());
    BOOST_CHECK_EQUAL(1, world.list_particles(species).size());
}

BOOST_AUTO_TEST_CASE(AddMolecule)
{
    const Species species("Species", /* radius = */ "1.0e-8", /* D = */ "1e-12");
    coordinate_type coordinate(50000);

    BOOST_CHECK(world.new_voxel(species, coordinate).second);
    BOOST_CHECK_EQUAL(1, world.num_particles(species));

    BOOST_CHECK(!world.get_voxel_pool_at(coordinate)->is_vacant());
}

BOOST_AUTO_TEST_CASE(AddMolecules)
{
    const Species species("Species", /* radius = */ "1.0e-8", /* D = */ "1e-12");
    const Integer N(60);

    BOOST_CHECK(world.add_molecules(species, N));
    BOOST_CHECK_EQUAL(N, world.num_particles(species));
}

BOOST_AUTO_TEST_CASE(Neighbor)
{
    const coordinate_type coordinate(26 + 52 * 26 + 52 * 52 * 26);
    const Real3 center(world.coordinate2position(coordinate));

    const Integer num_neighbors(world.num_neighbors(coordinate));
    for (Integer i(0); i < num_neighbors; ++i)
    {
        const coordinate_type neighbor(world.get_neighbor(coordinate, i));
        const Real3 position(world.coordinate2position(neighbor));
        BOOST_CHECK(length(position-center) < VOXEL_RADIUS*2.1);
    }
}

BOOST_AUTO_TEST_CASE(Move)
{
    const Species species("Species", /* radius = */ "1.0e-8", /* D = */ "1e-12");
    const coordinate_type src(world.inner2coordinate(10));
    const coordinate_type dest(world.inner2coordinate(1000));

    BOOST_CHECK(world.new_voxel(species, src).second);
    BOOST_CHECK(!world.get_voxel_pool_at(src)->is_vacant());
    BOOST_CHECK(world.get_voxel_pool_at(dest)->is_vacant());

    BOOST_CHECK(world.move(src, dest));
    BOOST_CHECK(world.get_voxel_pool_at(src)->is_vacant());
    BOOST_CHECK(!world.get_voxel_pool_at(dest)->is_vacant());
}

BOOST_AUTO_TEST_CASE(AddShape)
{
    const Species species("Species", /* radius = */ "1.0e-8", /* D = */ "1e-12");
    boost::shared_ptr<const Sphere> sphere(
            new Sphere(Real3(.5e-6, .5e-6, .5e-6), .45e-6));

    const Integer N(world.add_structure(species, sphere));
    BOOST_CHECK(N > 0);
    BOOST_CHECK_EQUAL(N, world.num_particles(species));
}

BOOST_AUTO_TEST_CASE(Structure)
{
    const Species
        membrane("Membrane", /* radius = */ "1.0e-8", /* D = */ "0"),
        species("Species", /* radius = */ "1.0e-8", /* D = */ "1e-12", /* location = */ "Membrane");
    boost::shared_ptr<const Sphere> sphere(
            new Sphere(Real3(.5e-6, .5e-6, .5e-6), .45e-6));

    BOOST_CHECK(world.add_structure(membrane, sphere) > 0);
    BOOST_CHECK(world.new_particle(
                Particle(species,
                         /* position = */ Real3(.5e-6, .5e-6, .95e-6),
                         /*   radius = */ 1e-8,
                         /*        D = */ 1e-12)).second);
}

BOOST_AUTO_TEST_SUITE_END()
