#include <iostream>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include <ecell4/spatiocyte/SpatiocyteSimulator.hpp>

using namespace ecell4;

typedef spatiocyte::SpatiocyteWorld world_type;
typedef spatiocyte::SpatiocyteSimulator simulator_type;
typedef VoxelSpaceBase::coordinate_type coordinate_type;

void
add_plane_offlattice_space(boost::shared_ptr<world_type> world,
                           const Real3& e1,
                           const Real3& e2,
                           const Species& species)
{
    OffLatticeSpace::position_container positions;
    OffLatticeSpace::coordinate_pair_list_type adjoining_pairs;

    const Real voxel_radius(world->voxel_radius());

    const Real len1(length(e1));
    const Real len2(length(e2));
    const int max1(len1 / voxel_radius);
    const int max2(len2 / voxel_radius);

    positions.reserve(max1*max2);

    for (int i(0); i < max1; ++i)
    {
        const Real3 v1(e1 * i * voxel_radius / len1);
        for (int j(0); j < max2; ++j)
            positions.push_back(v1 + e2 * j * voxel_radius / len2);
    }

    adjoining_pairs.reserve(max1*(max2-1) + (max1-1)*max2);

    for (int i(0); i < max1; ++i)
        for (int j(1); j < max2; ++j)
            adjoining_pairs.push_back(std::make_pair(i*max2+j-1, i*max2+j));

    for (int i(1); i < max1; ++i)
        for (int j(0); j < max2; ++j)
            adjoining_pairs.push_back(std::make_pair((i-1)*max2+j, i*max2+j));

    add_offlattice_space(world.get(), positions, adjoining_pairs, species);
}

int main(int argc, char** argv)
{
    const Real world_size(1e-6);
    const Real3 edge_lengths(world_size, world_size, world_size);
    const Real volume(world_size * world_size * world_size);
    const Real voxel_radius(2.5e-9);

    const Integer N(10000);

    const std::string radius("2.5e-9");

    const Real kd(1.0e+1);
    const Real ka(1.0e-4);

    Species speciesA("A", radius, "1e-9"),
            speciesB("B", radius, "0.0"),
            speciesC("C", radius, "1e-9");

    speciesC.set_attribute("location", "B");

    ReactionRule rr1(create_binding_reaction_rule(speciesA, speciesB, speciesC, ka));
    ReactionRule rr2(create_unbinding_reaction_rule(speciesC, speciesA, speciesB, kd));

    boost::shared_ptr<NetworkModel> model(new NetworkModel());
    model->add_species_attribute(speciesA);
    model->add_species_attribute(speciesB);
    model->add_species_attribute(speciesC);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    boost::shared_ptr<GSLRandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    rng->seed(1);

    boost::shared_ptr<world_type> world(new world_type(edge_lengths, voxel_radius, rng));

    add_plane_offlattice_space(
            world,
            Real3(world_size, 0, world_size/2),
            Real3(0, world_size, world_size/2),
            speciesB);
    world->add_molecules(speciesA, N);

    simulator_type sim(model, world);
    Real next_time(0.0), dt(0.00001);
    std::cout << "t"
              << "\t" << speciesA.serial()
              << "\t" << speciesB.serial()
              << "\t" << speciesC.serial()
              << std::endl;
    std::cout << sim.t()
              << "\t" << world->num_molecules(speciesA)
              << "\t" << world->num_molecules(speciesB)
              << "\t" << world->num_molecules(speciesC)
              << std::endl;
    for (int i(0); i < 100; ++i)
    {
        next_time += dt;
        while (sim.step(next_time)) {}

        std::cout << sim.t()
                  << "\t" << world->num_molecules(speciesA)
                  << "\t" << world->num_molecules(speciesB)
                  << "\t" << world->num_molecules(speciesC)
                  << std::endl;
    }

    return 0;
}
