#define BOOST_TEST_MODULE "SpatiocyteDiffusion_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/observers.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include "../SpatiocyteSimulator.hpp"

using namespace ecell4;

typedef spatiocyte::SpatiocyteWorld world_type;
typedef spatiocyte::SpatiocyteSimulator simulator_type;

BOOST_AUTO_TEST_CASE(Diffusion3D)
{
    const Real  VOXEL_RADIUS(2.5e-9);
    const Real3 EDGE_LENGTHS(0.5e-6, 0.5e-6, 0.5e-6);

    boost::shared_ptr<GSLRandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    boost::shared_ptr<NetworkModel> model(new NetworkModel());

    const Species species("A", /* radius = */ "2.5e-9", /* D = */ "1e-12");
    model->add_species_attribute(species);

    boost::shared_ptr<world_type> world(new world_type(
                new LatticeSpaceVectorImpl(EDGE_LENGTHS, VOXEL_RADIUS),
                rng));
    world->bind_to(model);
    world->add_molecules(species, 60);

    simulator_type simulator(model, world);

    boost::shared_ptr<FixedIntervalTrajectoryObserver> observer(
            new FixedIntervalTrajectoryObserver(/* dt = */ 1.0e-4));

    simulator.run(1.0, observer);

    Integer num_of_data(0);
    Real rss(0.0);

    for (Integer step(0); step < observer->count(); ++step)
    {
        Integer count(0);
        Real sum(0.0);

        for (Integer i(0); i < observer->num_tracers(); ++i)
        {
            const std::vector<Real3>& trajectory(observer->data().at(i));
            if (step < trajectory.size())
            {
                count++;
                sum += length_sq(trajectory.at(step) - trajectory.at(0));
            }
        }

        if (count != 0)
        {
            const Real msd(sum / observer->num_tracers());
            const Real correct_msd(6.0 * 1.0e-12 * observer->t().at(step));

            // std::cout << observer->t().at(step) << " " << msd << " " << correct_msd << std::endl;

            num_of_data++;
            rss += pow(msd - correct_msd, 2);
        }
    }

    BOOST_CHECK(rss / num_of_data < 1.0e-20);
}
