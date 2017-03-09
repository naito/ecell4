#define BOOST_TEST_MODULE "SpatiocyteSimulator_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Sphere.hpp>

#include <ecell4/core/LatticeSpaceVectorImpl.hpp>
#include <ecell4/core/LatticeSpaceCellListImpl.hpp>
#include <ecell4/core/OffLatticeSpace.hpp>

#include "../SpatiocyteSimulator.hpp"

using namespace ecell4;
using namespace ecell4::spatiocyte;

const Real DEFAULT_VOXEL_RADIUS = 1e-8;

struct VectorImplFixture
{
    const Real VOXEL_RADIUS;
    boost::shared_ptr<NetworkModel> model;
    boost::shared_ptr<GSLRandomNumberGenerator> rng;
    boost::shared_ptr<SpatiocyteWorld> world;
    SpatiocyteSimulator simulator;

    VectorImplFixture()
        : VOXEL_RADIUS(1e-8),
          model(new NetworkModel()),
          rng(new GSLRandomNumberGenerator()),
          world(new SpatiocyteWorld(
                      new LatticeSpaceVectorImpl( /* edge_lengths = */ Real3(1e-6, 1e-6, 1e-6),
                          /* voxel_radius = */ VOXEL_RADIUS),
                      rng)),
          simulator(model, world)
    {
        world->bind_to(model);
    }
};

BOOST_FIXTURE_TEST_SUITE(VectorImplTest, VectorImplFixture)


BOOST_AUTO_TEST_CASE(constructor) {}

BOOST_AUTO_TEST_CASE(StepWithSingleParticle)
{
    const ecell4::Species species("A", /* radius = */ "2.5e-9", /* D = */ "1e-12");
    model->add_species_attribute(species);

    BOOST_ASSERT(world->new_voxel(species, world->inner2coordinate(0)).second);
    BOOST_CHECK_EQUAL(1, world->num_voxels(species));

    simulator.initialize();
    for (int i(0); i < 10; ++i)
    {
        BOOST_CHECK_NO_THROW(simulator.step());
    }

    BOOST_CHECK_EQUAL(1, world->num_voxels(species));
}

BOOST_AUTO_TEST_CASE(StepWithSingleSpecies)
{
    const ecell4::Species species("A", /* radius = */ "2.5e-9", /* D = */ "1e-12");
    model->add_species_attribute(species);

    BOOST_ASSERT(world->add_molecules(species, 60));
    BOOST_CHECK_EQUAL(60, world->num_voxels(species));

    simulator.initialize();
    for (int i(0); i < 10; ++i)
    {
        BOOST_CHECK_NO_THROW(simulator.step());
    }

    BOOST_CHECK_EQUAL(60, world->num_voxels(species));
}

BOOST_AUTO_TEST_CASE(UnimolecularReaction)
{
    const ecell4::Species
        speciesA("A", /* radius = */ "2.5e-9", /* D = */ "1e-12"),
        speciesB("B", /* radius = */ "2.5e-9", /* D = */ "1e-12");

    model->add_species_attribute(speciesA);
    model->add_species_attribute(speciesB);

    model->add_reaction_rule(
            create_unimolecular_reaction_rule(
                /*    reactant = */ speciesA,
                /*     product = */ speciesB,
                /* coefficient = */ 1e6));

    BOOST_ASSERT(world->add_molecules(speciesA, 25));
    BOOST_CHECK_EQUAL(25, world->num_voxels(speciesA));
    BOOST_CHECK_EQUAL( 0, world->num_voxels(speciesB));

    simulator.initialize();
    for (int i(0); i < 10; ++i)
    {
        BOOST_CHECK_NO_THROW(simulator.step());
    }

    BOOST_CHECK(world->num_voxels(speciesB) > 0);
    BOOST_CHECK_EQUAL(25, world->num_voxels(speciesA) + world->num_voxels(speciesB));
}

BOOST_AUTO_TEST_CASE(BindingReaction)
{
    const ecell4::Species
        speciesA("A", /* radius = */ "2.5e-9", /* D = */ "1e-12"),
        speciesB("B", /* radius = */ "2.5e-9", /* D = */ "1e-12"),
        speciesC("C", /* radius = */ "2.5e-9", /* D = */ "1e-12");

    model->add_species_attribute(speciesA);
    model->add_species_attribute(speciesB);
    model->add_species_attribute(speciesC);

    model->add_reaction_rule(
            create_binding_reaction_rule(
                /* reactant[0] = */ speciesA,
                /* reactant[1] = */ speciesB,
                /*  product    = */ speciesC,
                /* coefficient = */ 1e-20));

    BOOST_ASSERT(world->add_molecules(speciesA, 25));
    BOOST_ASSERT(world->add_molecules(speciesB, 25));
    BOOST_CHECK_EQUAL(25, world->num_voxels(speciesA));
    BOOST_CHECK_EQUAL(25, world->num_voxels(speciesB));
    BOOST_CHECK_EQUAL( 0, world->num_voxels(speciesC));

    simulator.initialize();
    while (!simulator.check_reaction())
    {
        BOOST_CHECK_NO_THROW(simulator.step());
    }

    const Integer numberOfC(world->num_voxels(speciesC));
    BOOST_CHECK(numberOfC > 0);
    BOOST_CHECK_EQUAL(25, world->num_voxels(speciesA) + numberOfC);
    BOOST_CHECK_EQUAL(25, world->num_voxels(speciesB) + numberOfC);
}

BOOST_AUTO_TEST_CASE(UnbindingReaction)
{
    const ecell4::Species
        speciesA("A", /* radius = */ "2.5e-9", /* D = */ "1e-12"),
        speciesB("B", /* radius = */ "2.5e-9", /* D = */ "1e-12"),
        speciesC("C", /* radius = */ "2.5e-9", /* D = */ "1e-12");

    model->add_species_attribute(speciesA);
    model->add_species_attribute(speciesB);
    model->add_species_attribute(speciesC);

    model->add_reaction_rule(
            create_unbinding_reaction_rule(
                /* reactant    = */ speciesA,
                /*  product[0] = */ speciesB,
                /*  product[1] = */ speciesC,
                /* coefficient = */ 1e5));

    BOOST_ASSERT(world->add_molecules(speciesA, 25));
    BOOST_CHECK_EQUAL(25, world->num_voxels(speciesA));

    simulator.initialize();
    while (!simulator.check_reaction())
    {
        BOOST_CHECK_NO_THROW(simulator.step());
    }

    const Integer numberOfA(world->num_voxels(speciesA));
    BOOST_CHECK(numberOfA > 0);
    BOOST_CHECK_EQUAL(25, world->num_voxels(speciesB) + numberOfA);
    BOOST_CHECK_EQUAL(25, world->num_voxels(speciesC) + numberOfA);
}

BOOST_AUTO_TEST_CASE(DegradationReaction)
{
    const ecell4::Species species("A", /* radius = */ "2.5e-9", /* D = */ "1e-12");
    model->add_species_attribute(species);

    model->add_reaction_rule(
            create_degradation_reaction_rule(
                /*     species = */ species,
                /* coefficient = */ 1e5));

    BOOST_ASSERT(world->add_molecules(species, 25));

    simulator.initialize();
    while (!simulator.check_reaction())
    {
        BOOST_CHECK_NO_THROW(simulator.step());
    }

    BOOST_CHECK(world->num_voxels(species) < 25);
}

BOOST_AUTO_TEST_CASE(Scheduler)
{
    const ecell4::Species
        speciesA("A", /* radius = */ "2.5e-9", /* D = */ "1.0e-12"),
        speciesB("B", /* radius = */ "2.5e-9", /* D = */ "1.1e-12"),
        speciesC("C", /* radius = */ "2.5e-9", /* D = */ "1.2e-12");

    model->add_species_attribute(speciesA);
    model->add_species_attribute(speciesB);
    model->add_species_attribute(speciesC);

    BOOST_ASSERT(world->new_voxel(speciesA, world->inner2coordinate(0)).second);
    BOOST_ASSERT(world->new_voxel(speciesB, world->inner2coordinate(1)).second);
    BOOST_ASSERT(world->new_voxel(speciesC, world->inner2coordinate(2)).second);

    BOOST_ASSERT(world->has_molecule_pool(speciesA));
    BOOST_ASSERT(world->has_molecule_pool(speciesB));
    BOOST_ASSERT(world->has_molecule_pool(speciesC));

    const MoleculePool
        *mpA(world->find_molecule_pool(speciesA)),
        *mpB(world->find_molecule_pool(speciesB)),
        *mpC(world->find_molecule_pool(speciesC));

    BOOST_ASSERT(mpA->begin() != mpA->end());
    BOOST_ASSERT(mpB->begin() != mpB->end());
    BOOST_ASSERT(mpC->begin() != mpC->end());

    SpatiocyteWorld::coordinate_type
        coordinateA((*mpA->begin()).coordinate),
        coordinateB((*mpB->begin()).coordinate),
        coordinateC((*mpC->begin()).coordinate);

    simulator.initialize();

    simulator.step();
    BOOST_ASSERT((*mpA->begin()).coordinate == coordinateA);
    BOOST_ASSERT((*mpB->begin()).coordinate == coordinateB);
    BOOST_ASSERT((*mpC->begin()).coordinate != coordinateC);

    coordinateC = (*mpC->begin()).coordinate;

    simulator.step();
    BOOST_ASSERT((*mpA->begin()).coordinate == coordinateA);
    BOOST_ASSERT((*mpB->begin()).coordinate != coordinateB);
    BOOST_ASSERT((*mpC->begin()).coordinate == coordinateC);

    coordinateB = (*mpB->begin()).coordinate;

    simulator.step();
    BOOST_ASSERT((*mpA->begin()).coordinate != coordinateA);
    BOOST_ASSERT((*mpB->begin()).coordinate == coordinateB);
    BOOST_ASSERT((*mpC->begin()).coordinate == coordinateC);
}


BOOST_AUTO_TEST_SUITE_END()
