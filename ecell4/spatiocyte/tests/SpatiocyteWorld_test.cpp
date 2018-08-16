#define BOOST_TEST_MODULE "SpatiocyteWorld_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include "../SpatiocyteWorld.hpp"
#include "../../core/Sphere.hpp"
//#include <ecell4/core/Sphere.hpp>
#include <fstream>

using namespace ecell4;
using namespace ecell4::spatiocyte;

const Real DEFAULT_VOXEL_RADIUS = 1e-8;

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_constructor)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Real3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_t)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Real3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);
    BOOST_CHECK_EQUAL(world.t(), 0);
    world.set_t(23.4);
    BOOST_CHECK_EQUAL(world.t(), 23.4);
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_num_species)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Real3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);
    BOOST_CHECK_EQUAL(world.list_species().size(), 0);
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_has_species)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Real3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);
    Species sp(std::string("Species"));
    BOOST_CHECK(!world.has_species(sp));
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_list_particles)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    const Real3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);
    std::vector<std::pair<ParticleID, Particle> > particles(world.list_particles());
    BOOST_CHECK_EQUAL(particles.size(), 0);
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_update_particles)
{
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SerialIDGenerator<ParticleID> sidgen;
    const Real3 edge_lengths(1e-6, 1e-6, 1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);

    ParticleID pid(sidgen());
    Species sp(std::string("A"));
    const Real3 pos(2e-7, 1e-7, 0);
    Real r(0);
    Real d(0);
    Particle p(sp, pos, r, d);

    world.update_particle(pid, p);

    BOOST_CHECK(world.has_species(sp));
    BOOST_CHECK(world.has_particle(pid));
    BOOST_CHECK_EQUAL(world.list_particles().size(), 1);
    BOOST_CHECK_EQUAL(world.list_particles(sp).size(), 1);
}

// BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_register_species)
// {
//     const Real3 edge_lengths(1e-6,1e-6,1e-6);
//     const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
//     boost::shared_ptr<GSLRandomNumberGenerator>
//         rng(new GSLRandomNumberGenerator());
//     SpatiocyteWorld world(edge_lengths, voxel_radius, rng);
// 
//     Species sp(std::string("TEST"));
// 
//     BOOST_CHECK(world.register_species(sp));
//     BOOST_CHECK(world.has_species(sp));
// 
//     std::vector<Species> list;
//     list.push_back(sp);
// 
//     BOOST_CHECK(list == world.list_species());
// }

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_add_molecule)
{
    const Real3 edge_lengths(1e-6,1e-6,1e-6);
    const Real voxel_radius(2.5e-9);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");

    const Voxel voxel(world.position2voxel(edge_lengths / 2.0));
    // BOOST_CHECK(world.place_voxel(sp, coord).second);
    BOOST_CHECK(world.new_voxel(sp, voxel));
    BOOST_CHECK_EQUAL(world.num_particles(sp), 1);

    const VoxelPool* mt(voxel.get_voxel_pool());
    BOOST_CHECK(!mt->is_vacant());
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_add_molecules)
{
    const Real3 edge_lengths(1e-6,1e-6,1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");
    const Integer N(60);

    BOOST_CHECK(world.add_molecules(sp, N));
    BOOST_CHECK_EQUAL(world.num_particles(sp), N);
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_neighbor)
{
    const Real3 edge_lengths(1e-6,1e-6,1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);

    const Voxel voxel(world.position2voxel(edge_lengths / 2.0));
    const Real3 cp(voxel.position());

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");
    for (Integer i(0); i < voxel.num_neighbors(); ++i)
    {
        world.new_voxel(sp, voxel.get_neighbor(i));
    }
    std::vector<std::pair<ParticleID, Particle> > particles(world.list_particles());
    for (std::vector<std::pair<ParticleID, Particle> >::iterator itr(
                particles.begin()); itr != particles.end(); ++itr)
    {
        Real3 pos((*itr).second.position());
        BOOST_ASSERT(length(pos-cp) < voxel_radius*2.1);
    }

#ifdef WITH_HDF5
    world.save("neighbor.h5");
#endif
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_add_shape)
{
    const Real3 edge_lengths(1e-6,1e-6,1e-6);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");

    boost::shared_ptr<const Sphere> sphere(new Sphere(Real3(5e-7, 5e-7, 5e-7), 5e-7*1.5));

    const Integer n(world.add_structure(sp, sphere));
    BOOST_ASSERT(n > 0);
    BOOST_CHECK_EQUAL(world.num_particles(sp), n);

#ifdef WITH_HDF5
    world.save("sphere.h5");
#endif
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_move)
{
    const Real3 edge_lengths(1e-6,1e-6,1e-6);
    const Real voxel_radius(2.5e-9);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);

    Species sp(std::string("TEST"));
    sp.set_attribute("radius", "2.5e-9");
    sp.set_attribute("D", "1e-12");

    const Voxel from(world.position2voxel(Real3(0.3e-6, 0.5e-6, 0.5e-6)));
    const Voxel to(world.position2voxel(Real3(0.5e-6, 0.5e-6, 0.5e-6)));

    BOOST_CHECK(world.new_voxel(sp, from));
    BOOST_CHECK(world.move(from, to));

    const VoxelPool* mt(to.get_voxel_pool());
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(world.move(from, to));
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_structure)
{
    const Real3 edge_lengths(5e-7, 5e-7, 5e-7);
    const Real voxel_radius(DEFAULT_VOXEL_RADIUS);
    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    SpatiocyteWorld world(edge_lengths, voxel_radius, rng);

    Species membrane("Membrane", "2.5e-9", "0");

    Species sp("SpeciesA", "2.5e-9", "1e-12");
    sp.set_attribute("location", "Membrane");

    boost::shared_ptr<const Sphere> sphere(new Sphere(Real3(2.5e-7, 2.5e-7, 2.5e-7), 2e-7));

    BOOST_CHECK(world.add_structure(membrane, sphere) == 5892);
    BOOST_CHECK(!world.new_particle(Particle(sp, Real3(2.5e-7, 2.5e-7, 4.5e-7), 2.5e-9, 1e-12)));
    BOOST_CHECK(world.new_particle(Particle(sp, Real3(2.5e-7, 2.5e-7, 4.5e-7 - voxel_radius * 2), 2.5e-9, 1e-12)));

#ifdef WITH_HDF5
    world.save("structure.h5");
#endif
}
