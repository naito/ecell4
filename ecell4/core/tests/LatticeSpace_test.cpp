#define BOOST_TEST_MODULE "LatticeSpace_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/MolecularType.hpp>
#include <ecell4/core/LatticeSpaceVectorImpl.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>

using namespace ecell4;

typedef VoxelSpaceBase::coordinate_type coordinate_type;

struct LatticeSpaceVectorFixture
{
    const Real3 edge_lengths;
    const Real voxel_radius;
    LatticeSpaceVectorImpl space;
    SerialIDGenerator<ParticleID> sidgen;
    const Real D, radius;
    const Species sp;
    LatticeSpaceVectorFixture() :
        edge_lengths(2.5e-8, 2.5e-8, 2.5e-8),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, false),
        sidgen(), D(1e-12), radius(2.5e-9),
        sp("A", "2.5e-9", "1e-12")
    {}

    inline
    const Voxel gen_voxel(const Species& species, const Real& x, const Real& y, const Real& z)
    {
        return Voxel(species, space.position2coordinate(Real3(x, y, z)), radius, D);
    }
};

BOOST_FIXTURE_TEST_SUITE( LatticeSpaceVectorTest, LatticeSpaceVectorFixture )

BOOST_AUTO_TEST_CASE( TestConstructor ) {}

BOOST_AUTO_TEST_CASE( TestValues )
{
    BOOST_CHECK( space.inner_size() < space.size() );
}

BOOST_AUTO_TEST_CASE( TestVoxel )
{
    BOOST_CHECK_EQUAL( space.num_voxels(),  0 );
    BOOST_CHECK_EQUAL( space.num_voxels(sp), 0 );
    BOOST_CHECK_EQUAL( space.num_voxels_exact(sp), 0 );

    BOOST_CHECK( ! space.has_voxel(ParticleID()) );
    BOOST_CHECK( ! space.has_voxel(ParticleID(std::make_pair(0, 1))));

    BOOST_CHECK_THROW( space.get_voxel(ParticleID()), NotFound );
    BOOST_CHECK_THROW( space.get_voxel(ParticleID(std::make_pair(0, 1))), NotFound );

    BOOST_CHECK_EQUAL( space.list_voxels().size(), 0);
    BOOST_CHECK_EQUAL( space.list_voxels(sp).size(), 0);
    BOOST_CHECK_EQUAL( space.list_voxels_exact(sp).size(), 0);
}

BOOST_AUTO_TEST_CASE( TestVoxelPool )
{
    BOOST_CHECK( ! space.has_voxel_pool(sp) );
    BOOST_CHECK_THROW( space.find_voxel_pool(sp), NotFound);

    BOOST_CHECK( ! space.has_molecule_pool(sp) );
    BOOST_CHECK_THROW( space.find_molecule_pool(sp), NotFound);
}

BOOST_AUTO_TEST_CASE( TestSpecies )
{
    BOOST_CHECK_EQUAL( space.num_species(), 0 );
    BOOST_CHECK( ! space.has_species(sp) );
}

BOOST_AUTO_TEST_CASE( TestUpdateVoxel )
{
    const ParticleID pidA(sidgen());
    const ParticleID pidB(sidgen());
    const Species spB("SpeciesB");

    BOOST_CHECK( space.update_voxel(pidA, gen_voxel(sp,  2.0e-8, 1.7e-8, 1.5e-8)) );
    BOOST_CHECK( space.update_voxel(pidB, gen_voxel(spB, 1.0e-8, 2.0e-8, 0.1e-8)) );

    BOOST_CHECK( space.has_voxel(pidA) );
    BOOST_CHECK_NO_THROW( space.get_voxel(pidA) );
    BOOST_CHECK_EQUAL( space.list_voxels().size(), 2 );
    BOOST_CHECK_EQUAL( space.list_voxels(sp).size(), 1 );
    BOOST_CHECK_EQUAL( space.list_voxels(spB).size(), 1 );
    BOOST_CHECK_EQUAL( space.list_voxels_exact(sp).size(), 1 );
    BOOST_CHECK_EQUAL( space.list_voxels_exact(spB).size(), 1 );

    /* ### Issue ###
     *
     * has_voxel_pool() returns whether a given species is hold as a structure.
     * However find_voxel_pool() returns a structure pool or a molecule pool
     * corresponding to a given species.
     */
    BOOST_CHECK( ! space.has_voxel_pool(sp) );
    BOOST_CHECK_NO_THROW( space.find_voxel_pool(sp) );
    BOOST_CHECK( ! space.has_voxel_pool(spB) );
    BOOST_CHECK_NO_THROW( space.find_voxel_pool(spB) );

    BOOST_CHECK( space.has_molecule_pool(sp) );
    BOOST_CHECK_NO_THROW( space.find_voxel_pool(sp) );
    BOOST_CHECK( space.has_molecule_pool(spB) );
    BOOST_CHECK_NO_THROW( space.find_voxel_pool(spB) );

    BOOST_CHECK_EQUAL( space.num_species(), 2 );
    BOOST_CHECK( space.has_species(sp) );
    BOOST_CHECK( space.has_species(spB) );

    BOOST_CHECK_EQUAL( space.list_particles().size(), 2 );
    BOOST_CHECK_EQUAL( space.list_particles(sp).size(), 1 );
    BOOST_CHECK_EQUAL( space.list_particles(spB).size(), 1 );
}

BOOST_AUTO_TEST_CASE( TestCoordinateTranslation )
{
    const Integer size((space.col_size()+2) * (space.layer_size() + 2) * (space.row_size() + 2));

    for (Integer coord(0); coord < size; ++coord)
    {
        BOOST_CHECK_EQUAL( space.global2coordinate(space.coordinate2global(coord)), coord );
        BOOST_CHECK_EQUAL( space.position2coordinate(space.coordinate2position(coord)), coord );
    }

    const Real3 zero_pos(space.coordinate2position(0));
    BOOST_ASSERT( zero_pos[0] < 0 );
    BOOST_ASSERT( zero_pos[1] < 0 );
    BOOST_ASSERT( zero_pos[2] < 0 );

    const Real3 origin_pos(
            space.coordinate2position((space.col_size() + 3) * (space.row_size() + 2) + 1));
    BOOST_CHECK_EQUAL( origin_pos[0], 0 );
    BOOST_CHECK_EQUAL( origin_pos[1], 0 );
    BOOST_CHECK_EQUAL( origin_pos[2], 0 );
}

BOOST_AUTO_TEST_CASE( TestAddRemove )
{
    const coordinate_type coord(space.position2coordinate(Real3(2.0e-8, 1.5e-8, 0.2e-8)));

    BOOST_CHECK( space.update_voxel(sidgen(), gen_voxel(sp, 2.0e-8, 1.5e-8, 0.2e-8)) );
    BOOST_CHECK_EQUAL( space.num_voxels(sp), 1 );

    BOOST_CHECK( ! space.get_voxel_pool_at(coord)->is_vacant());

    BOOST_CHECK( space.remove_voxel(coord) );
    BOOST_CHECK( space.get_voxel_pool_at(coord)->is_vacant() );
}

BOOST_AUTO_TEST_CASE( TestMove )
{
    const coordinate_type src(space.global2coordinate(Integer3(3, 4, 5)));
    const coordinate_type dest(space.global2coordinate(Integer3(3, 5, 5)));

    BOOST_CHECK( space.update_voxel(sidgen(), Voxel(sp, src, radius, D)) );

    BOOST_CHECK( ! space.get_voxel_pool_at(src)->is_vacant() );
    BOOST_CHECK( space.get_voxel_pool_at(dest)->is_vacant() );

    BOOST_CHECK( space.move(src, dest) );

    BOOST_CHECK( space.get_voxel_pool_at(src)->is_vacant() );
    BOOST_CHECK( ! space.get_voxel_pool_at(dest)->is_vacant() );

    BOOST_CHECK( space.update_voxel(sidgen(), Voxel(sp, src, radius, D)) );
    BOOST_CHECK( ! space.move(src, dest) );
}

BOOST_AUTO_TEST_CASE( TestNeighbor )
{
    for (Integer inner_coord(0); inner_coord < space.inner_size(); ++inner_coord)
    {
        const coordinate_type coord(space.inner2coordinate(inner_coord));
        const Real3 center(space.coordinate2position(coord));

        BOOST_ASSERT( space.num_neighbors(coord) == 12 );

        for (int i(0); i < 12; ++i)
        {
            const coordinate_type neighbor(space.get_neighbor(coord, i));
            if (!space.is_inside(neighbor))
                continue;
            const Real3 pos(space.coordinate2position(neighbor));
            const Real r_ratio(length(pos-center)/voxel_radius/2);
            BOOST_CHECK( r_ratio < 1.0001 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()


struct PeriodicFixture
{
    const Real3 edge_lengths;
    const Real voxel_radius;
    LatticeSpaceVectorImpl space;
    SerialIDGenerator<ParticleID> sidgen;
    const Real D, radius;
    const Species sp;
    PeriodicFixture() :
        edge_lengths(2.5e-8, 2.5e-8, 2.5e-8),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, true),
        sidgen(), D(1e-12), radius(2.5e-9),
        sp(std::string("A"), "2.5e-9", "1e-12")
    {
    }
};

BOOST_FIXTURE_TEST_SUITE( PeriodicSpaceTest, PeriodicFixture )

BOOST_AUTO_TEST_CASE( TestCol )
{
    const Integer col_size(space.col_size());
    const Integer row_size(space.row_size());
    const Integer layer_size(space.layer_size());

    for (Integer i(0); i < row_size; ++i)
        for (Integer j(0); j < layer_size; ++j)
        {
            const coordinate_type coord(space.global2coordinate(Integer3(0, i, j)));

            BOOST_CHECK( space.update_voxel(sidgen(), Voxel(sp, coord, radius, D)) );
        }

    // from 0 to col_size-1
    for (Integer i(0); i < row_size; ++i)
        for (Integer j(0); j < layer_size; ++j)
        {
            const coordinate_type coord(space.global2coordinate(Integer3(0, i, j)));
            const Integer nrnd((j&1)==1?2:3);
            const coordinate_type neighbor(space.get_neighbor_boundary(coord, nrnd));

            BOOST_CHECK_EQUAL( space.coordinate2global(neighbor).col, col_size-1 );
            BOOST_CHECK( space.move(coord, neighbor) );
        }

    // from col_size-1 to 0
    for (Integer i(0); i < row_size; ++i)
        for (Integer j(0); j < layer_size; ++j)
        {
            const coordinate_type coord(space.global2coordinate(Integer3(col_size-1, i, j)));
            const Integer nrnd((j&1)==1?4:5);
            const coordinate_type neighbor(space.get_neighbor_boundary(coord, nrnd));

            BOOST_CHECK_EQUAL( space.coordinate2global(neighbor).col, 0 );
            BOOST_CHECK( space.move(coord, neighbor) );
        }
}

BOOST_AUTO_TEST_CASE( TestRow )
{
    const Integer col_size(space.col_size());
    const Integer row_size(space.row_size());
    const Integer layer_size(space.layer_size());

    for (Integer layer(0); layer < layer_size; ++layer)
        for (Integer col(0); col < col_size; ++col)
        {
            const coordinate_type coord(space.global2coordinate(Integer3(col, 0, layer)));

            BOOST_CHECK( space.update_voxel(sidgen(), Voxel(sp, coord, radius, D)) );
        }
    // from 0 to row_size-1
    int row(0);
    for (Integer layer(0); layer < layer_size; ++layer)
        for (Integer col(0); col < col_size; ++col)
        {
            const coordinate_type coord(space.global2coordinate(Integer3(col, row, layer)));
            const coordinate_type neighbor(space.get_neighbor_boundary(coord, /* rnd = */ 0));

            BOOST_CHECK_EQUAL( space.coordinate2global(neighbor).row, row_size-1 );
            BOOST_CHECK( space.move(coord, neighbor) );
        }
    // from row_size-1 to 0
    row = row_size - 1;
    for (Integer layer(0); layer < layer_size; ++layer)
        for (Integer col(0); col < col_size; ++col)
        {
            const coordinate_type coord(space.global2coordinate(Integer3(col, row, layer)));
            const coordinate_type neighbor(space.get_neighbor_boundary(coord, /* rnd = */ 1));

            BOOST_CHECK_EQUAL( space.coordinate2global(neighbor).row, 0 );
            BOOST_CHECK( space.move(coord, neighbor) );
        }
}

BOOST_AUTO_TEST_CASE( TestLayer )
{
    const Integer col_size(space.col_size());
    const Integer row_size(space.row_size());
    const Integer layer_size(space.layer_size());

    Integer layer(0);
    for (Integer row(0); row < row_size; ++row)
        for (Integer col(0); col < col_size; ++col)
        {
            const coordinate_type coord(space.global2coordinate(Integer3(col, row, layer)));

            BOOST_CHECK( space.update_voxel(sidgen(), Voxel(sp, coord, radius, D)) );
        }
    // from 0 to layer_size-1
    for (Integer row(0); row < row_size; ++row)
        for (Integer col(0); col < col_size; ++col)
        {
            const coordinate_type coord(space.global2coordinate(Integer3(col, row, layer)));
            const Integer nrnd((col&1)==1?8:9);
            const coordinate_type neighbor(space.get_neighbor_boundary(coord, nrnd));

            BOOST_CHECK_EQUAL( space.coordinate2global(neighbor).layer, layer_size-1 );
            BOOST_CHECK( space.move(coord, neighbor) );
        }
    // from layer_size-1 to 0
    layer = layer_size - 1;
    for (Integer row(0); row < row_size; ++row)
        for (Integer col(0); col < col_size; ++col)
        {
            const coordinate_type coord(space.global2coordinate(Integer3(col, row, layer)));
            const Integer nrnd((col&1)==1?10:11);
            const coordinate_type neighbor(space.get_neighbor_boundary(coord, nrnd));

            BOOST_CHECK_EQUAL( space.coordinate2global(neighbor).layer, 0 );
            BOOST_CHECK( space.move(coord, neighbor) );
        }
}

BOOST_AUTO_TEST_SUITE_END()


struct StructureFixture
{
    const Real3 edge_lengths;
    const Real voxel_radius;
    LatticeSpaceVectorImpl space;
    SerialIDGenerator<ParticleID> sidgen;
    const Real D, radius;
    const Species structure, sp;
    StructureFixture() :
        edge_lengths(2.5e-8, 2.5e-8, 2.5e-8),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, false),
        sidgen(), D(1e-12), radius(2.5e-9),
        structure("Structure", "2.5e-9", "0"),
        sp("A", "2.5e-9", "1e-12", "Structure")
    {
    }
};

BOOST_FIXTURE_TEST_SUITE( StructureTest, StructureFixture )

BOOST_AUTO_TEST_CASE( TestStructureUpdate )
{
    const Real3 pos(2.7e-9, 1.3e-8, 2.0e-8);

    BOOST_CHECK( space.update_structure(Particle(structure, pos, radius, D)) );
    BOOST_CHECK_EQUAL( space.list_particles().size(), 1 );

    ParticleID pid(sidgen());

    BOOST_CHECK( space.update_voxel(pid, Voxel(sp, space.position2coordinate(pos), radius, D, structure.serial())) );
    BOOST_CHECK_EQUAL( space.list_particles().size(),   1 );
    BOOST_CHECK_EQUAL( space.list_particles(sp).size(), 1 );

    BOOST_CHECK( space.remove_particle(pid) );
    BOOST_CHECK_EQUAL( space.list_particles().size(),   1 ); // TODO -> 0
    BOOST_CHECK_EQUAL( space.list_particles(sp).size(), 0 );

    Species sp2("B", "2.5e-9", "1e-12");
    BOOST_CHECK_THROW(
        space.update_voxel(sidgen(), Voxel(sp2, space.position2coordinate(pos), radius, D)),
        NotSupported
        );
}

BOOST_AUTO_TEST_CASE( TestStructureMove )
{
    const Real3 pos1(2.7e-9, 1.3e-8, 2.0e-8);
    const Real3 pos2(1.2e-8, 1.5e-8, 1.8e-8);

    BOOST_CHECK( space.update_structure(Particle(structure, pos1, radius, D)) );
    BOOST_CHECK_EQUAL( space.list_particles().size(), 1 );
    BOOST_CHECK( space.update_structure(Particle(structure, pos2, radius, D)) );
    BOOST_CHECK_EQUAL( space.list_particles().size(), 2 ); // TODO -> 0

    BOOST_CHECK( space.update_voxel(sidgen(), Voxel(sp, space.position2coordinate(pos1), radius, D, structure.serial())) );
    BOOST_CHECK_EQUAL( space.list_particles(sp).size(), 1 );
    BOOST_CHECK_EQUAL( space.list_particles(structure).size(), 1 );
    BOOST_CHECK_EQUAL( space.list_particles().size(), 2 ); // TODO -> 1

    const coordinate_type coord1(space.position2coordinate(pos1));
    const coordinate_type coord2(space.position2coordinate(pos2));

    BOOST_CHECK( space.move(coord1, coord2) );
    BOOST_CHECK_EQUAL( space.list_particles(sp).size(), 1 );
    BOOST_CHECK_EQUAL( space.list_particles(structure).size(), 1 );
    BOOST_CHECK_EQUAL( space.list_particles().size(), 2 ); // TODO -> 1
}

#ifdef WITH_HDF5
BOOST_AUTO_TEST_CASE( TestSaveAndLoad )
{
    space.make_structure_type(structure, Shape::TWO, "");
    const Integer layer(space.layer_size()/2);
    for (int col(0); col < space.col_size(); ++col)
    {
        for (int row(0); row < space.row_size(); ++row)
        {
            const Real3 pos(space.global2position(Integer3(col, row, layer)));
            BOOST_ASSERT( space.update_structure(Particle(structure, pos, radius, D)) );
        }
    }

    const coordinate_type center(space.global2coordinate(Integer3(space.col_size()/2, space.row_size()/2, layer)));

    BOOST_ASSERT( space.update_voxel(sidgen(), Voxel(sp, center, radius, D, structure.serial())) );
    // #XXX !!!Warning!!! Ideally, not necessary to give structure.serial() explicitly

    H5::H5File fout("data.h5", H5F_ACC_TRUNC);
    boost::scoped_ptr<H5::Group> group(new H5::Group(fout.createGroup("LatticeSpace")));
    space.save_hdf5(group.get());
    fout.close();

    LatticeSpaceVectorImpl space2(Real3(3e-8, 3e-8, 3e-8), voxel_radius);
    H5::H5File fin("data.h5", H5F_ACC_RDONLY);
    const H5::Group groupin(fin.openGroup("LatticeSpace"));
    space2.load_hdf5(groupin);
    fin.close();

    BOOST_CHECK_EQUAL( space.edge_lengths(),  space2.edge_lengths()  );
    BOOST_CHECK_EQUAL( space.voxel_radius(),  space2.voxel_radius()  );
    BOOST_CHECK_EQUAL( space.is_periodic(),   space2.is_periodic()   );
    BOOST_CHECK_EQUAL( space.t(),             space2.t()             );
    BOOST_CHECK_EQUAL( space.num_particles(), space2.num_particles() );
    BOOST_CHECK_EQUAL( space.num_species(),   space2.num_species()   );

    std::vector<Species> species(space.list_species());
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        const Species species((*itr).serial());

        const VoxelPool *vp1(space.find_voxel_pool(species));
        const VoxelPool *vp2(space2.find_voxel_pool(species));

        BOOST_CHECK_EQUAL( vp1->radius(),        vp2->radius()        );
        BOOST_CHECK_EQUAL( vp1->D(),             vp2->D()             );
        BOOST_CHECK_EQUAL( vp1->get_dimension(), vp2->get_dimension() );

        const MolecularType* mtb1(dynamic_cast<const MolecularType*>(vp1));
        const MolecularType* mtb2(dynamic_cast<const MolecularType*>(vp2));

        BOOST_ASSERT( (mtb1 && mtb2) || (!mtb1 && !mtb2) );

        if (!mtb1 || !mtb2)
        {
            continue;
        }

        MoleculePool::container_type voxels1, voxels2;
        std::copy(mtb1->begin(), mtb1->end(), back_inserter(voxels1));
        std::copy(mtb2->begin(), mtb2->end(), back_inserter(voxels2));

        BOOST_ASSERT( voxels1.size() == voxels2.size() );

        std::sort(voxels1.begin(), voxels1.end());
        std::sort(voxels2.begin(), voxels2.end());
        for (int i(0); i < voxels1.size(); ++i)
        {
            BOOST_CHECK_EQUAL( voxels1.at(i).pid,        voxels2.at(i).pid        );
            BOOST_CHECK_EQUAL( voxels1.at(i).coordinate, voxels2.at(i).coordinate );
        }
    }
}
#endif

BOOST_AUTO_TEST_SUITE_END()
