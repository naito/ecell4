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

OffLatticeSpace*
create_plane_offlattice(const Real& voxel_radius, const Real3& e1, const Real3& e2,
                        const Species& species) {

    OffLatticeSpace::position_container positions;
    OffLatticeSpace::coordinate_pair_list_type adjoining_pairs;

    const Real len1(length(e1));
    const Real len2(length(e2));
    const unsigned int max1(len1 / voxel_radius);
    const unsigned int max2(len2 / voxel_radius);

    positions.reserve(max1*max2);

    for (unsigned int i(0); i < max1; ++i)
    {
        const Real3 v1(e1 * i * voxel_radius / len1);
        for (unsigned int j(0); j < max2; ++j)
            positions.push_back(v1 + e2 * j * voxel_radius / len2);
    }

    adjoining_pairs.reserve(max1*(max2-1) + (max1-1)*max2);

    for (unsigned int i(0); i < max1; ++i)
        for (unsigned int j(1); j < max2; ++j)
            adjoining_pairs.push_back(std::make_pair(i*max2+j-1, i*max2+j));

    for (unsigned int i(1); i < max1; ++i)
        for (unsigned int j(0); j < max2; ++j)
            adjoining_pairs.push_back(std::make_pair((i-1)*max2+j, i*max2+j));

    OffLatticeSpace* space(new OffLatticeSpace(voxel_radius));
    space->reset(positions, adjoining_pairs);

    space->make_structure_type(species, Shape::TWO, "");
    const Real radius(species.get_attribute_as<Real>("radius"));

    for (unsigned int coord(0); coord < max1*max2; ++coord)
        space->update_voxel(ParticleID(), Voxel(species, coord, radius, 0.0, ""));

    return space;
}

int main(int argc, char** argv)
{
    const Real world_size(1e-6);
    const Real3 edge_lengths(world_size, world_size, world_size);
    const Real volume(world_size * world_size * world_size); const Real voxel_radius(2.5e-9);

    const Integer N(600);

    const std::string D("1e-12"), radius("2.5e-9");

    const Real kd(1.0e+1);
    const Real ka(1.0e-4);

    Species speciesA("A", radius, D),
            speciesB("B", radius, "0.0"),
            speciesC("C", radius, D);

    speciesC.set_attribute("location", "B");

    ReactionRule rr1(create_binding_reaction_rule(speciesA, speciesB, speciesC, ka));
    ReactionRule rr2(create_unbinding_reaction_rule(speciesC, speciesA, speciesB, kd));

    boost::shared_ptr<NetworkModel> model(new NetworkModel);
    model->add_species_attribute(speciesA);
    model->add_species_attribute(speciesB);
    model->add_species_attribute(speciesC);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    boost::shared_ptr<GSLRandomNumberGenerator> rng(new GSLRandomNumberGenerator());
    rng->seed(1);

    boost::shared_ptr<world_type> world(new world_type(edge_lengths, voxel_radius, rng));

    world->add_space(create_plane_offlattice(voxel_radius,
                Real3(world_size, 0, world_size/2),
                Real3(0, world_size, world_size/2),
                speciesB));
    world->add_molecules(speciesA, N);

    simulator_type sim(model, world);
    Real next_time(0.0), dt(0.02);
    std::cout << sim.t()
              << "\t" << world->num_molecules(speciesA)
              << "\t" << world->num_molecules(speciesB)
              << "\t" << world->num_molecules(speciesC)
              << std::endl;
    for (unsigned int i(0); i < 100; ++i)
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
