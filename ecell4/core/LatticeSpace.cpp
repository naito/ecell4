#include "LatticeSpace.hpp"
#include "hcp_lattice.hpp"

#ifdef WIN32_MSC
double rint(const double x);
double round(const double x);
#endif

namespace ecell4
{

LatticeSpace::LatticeSpace(const Real3& edge_lengths,
                           const Real& voxel_radius,
                           const bool is_periodic)
    : base_type(voxel_radius),
      edge_lengths_(edge_lengths)
{
    set_lattice_properties(is_periodic);
}

LatticeSpace::~LatticeSpace() {}

void LatticeSpace::set_lattice_properties(const bool is_periodic)
{
    HCP_L = voxel_radius_ / sqrt(3.0);
    HCP_X = voxel_radius_ * sqrt(8.0 / 3.0); // Lx
    HCP_Y = voxel_radius_ * sqrt(3.0); // Ly

    const Real lengthX = edge_lengths_[0];
    const Real lengthY = edge_lengths_[1];
    const Real lengthZ = edge_lengths_[2];

    col_size_   = (Integer)rint(lengthX / HCP_X) + 1;
    layer_size_ = (Integer)rint(lengthY / HCP_Y) + 1;
    row_size_   = (Integer)rint((lengthZ / 2) / voxel_radius_) + 1;

    if (is_periodic)
    {
        // The number of voxels in each axis must be even for a periodic boundary.
        col_size_ = (col_size_ % 2 == 0 ? col_size_ : col_size_ + 1);
        layer_size_ = (layer_size_ % 2 == 0 ? layer_size_ : layer_size_ + 1);
        row_size_ = (row_size_ % 2 == 0 ? row_size_ : row_size_ + 1);
    }

    row_size_ += 2;
    layer_size_ += 2;
    col_size_ += 2;
}

void LatticeSpace::reset(const Real3& edge_lengths,
                         const Real& voxel_radius,
                         const bool is_periodic)
{
    edge_lengths_ = edge_lengths;
    voxel_radius_ = voxel_radius;

    set_lattice_properties(is_periodic);
}

const Integer LatticeSpace::col_size() const
{
    return col_size_ - 2;
}

const Integer LatticeSpace::row_size() const
{
    return row_size_ - 2;
}

const Integer LatticeSpace::layer_size() const
{
    return layer_size_ - 2;
}

Integer3 LatticeSpace::coordinate2global(const coordinate_type& coord) const
{
    const Integer NUM_COLROW(row_size_ * col_size_);
    const Integer LAYER(coord / NUM_COLROW);
    const Integer SURPLUS(coord - LAYER * NUM_COLROW);
    const Integer COL(SURPLUS / row_size_);
    const Integer3 global(COL, SURPLUS - COL * row_size_, LAYER);
    const Integer3 retval(
        global.col - 1, global.row - 1, global.layer - 1);
    return retval;
}

LatticeSpace::coordinate_type
LatticeSpace::global2coordinate(const Integer3& global) const
{
    const Integer3 g(global.col + 1, global.row + 1, global.layer + 1);
    return g.row + row_size_ * (g.col + col_size_ * g.layer);
}

Real3 LatticeSpace::global2position(const Integer3& global) const
{
    return calc_lattice_position(voxel_radius_, global);
}

Integer3 LatticeSpace::position2global(const Real3& pos) const
{
    const Integer col(round(pos[0] / HCP_X));
    const Integer layer(round((pos[1] - (col % 2) * HCP_L) / HCP_Y));
    const Integer row(round((pos[2] / voxel_radius_ - ((layer + col) % 2)) / 2));
    const Integer3 global(col, row, layer);
    return global;
}

LatticeSpace::coordinate_type
LatticeSpace::periodic_transpose(const coordinate_type& coord) const
{
    Integer3 global(coordinate2global(coord));

    global.col = global.col % col_size();
    global.row = global.row % row_size();
    global.layer = global.layer % layer_size();

    global.col = global.col < 0 ? global.col + col_size() : global.col;
    global.row = global.row < 0 ? global.row + row_size() : global.row;
    global.layer = global.layer < 0 ? global.layer + layer_size() : global.layer;

    return global2coordinate(global);
}

bool LatticeSpace::is_in_range(const coordinate_type& coord) const
{
    return coord >= 0 && coord < row_size_ * col_size_ * layer_size_;
}

bool LatticeSpace::is_inside(const coordinate_type& coord) const
{
    const Integer3 global(coordinate2global(coord));
    return 0 <= global.col   && global.col   < col_size()
        && 0 <= global.row   && global.row   < row_size()
        && 0 <= global.layer && global.layer < layer_size();
}

/*
 * for LatticeSpaceBase
 */
Real3 LatticeSpace::coordinate2position(const coordinate_type& coord) const
{
    return global2position(coordinate2global(coord));
}

LatticeSpace::coordinate_type
LatticeSpace::position2coordinate(const Real3& pos) const
{
    return global2coordinate(position2global(pos));
}

Integer LatticeSpace::num_neighbors(const coordinate_type& coord) const
{
    if (!is_inside(coord)) return 0;
    return 12;
}

LatticeSpace::coordinate_type
LatticeSpace::get_neighbor(const coordinate_type& coord,
                           const Integer& nrand) const
{
    if (!is_inside(coord))
        throw NotFound("There is no neighbor voxel.");

    return calc_neighbor(Integer3(col_size_, row_size_, layer_size_), coord, nrand);
}

Integer LatticeSpace::size() const
{
    return row_size_ * col_size_ * layer_size_;
}

Integer3 LatticeSpace::shape() const
{
    return Integer3(col_size_, row_size_, layer_size_);
}

LatticeSpace::coordinate_type
LatticeSpace::inner2coordinate(const coordinate_type inner) const {
    const Integer num_row(row_size());
    const Integer num_col(col_size());

    const Integer NUM_COLROW(num_row * num_col);
    const Integer LAYER(inner / NUM_COLROW);
    const Integer SURPLUS(inner - LAYER * NUM_COLROW);
    const Integer COL(SURPLUS / num_row);
    const Integer3 g(COL, SURPLUS - COL * num_row, LAYER);

    return global2coordinate(g);
}

Integer LatticeSpace::inner_size() const
{
    return col_size() * row_size() * layer_size();
}

Real LatticeSpace::unit_voxel_volume() const
{
    return 4.0 * sqrt(2.0);
}

const Real3& LatticeSpace::edge_lengths() const
{
    return edge_lengths_;
}

Real3 LatticeSpace::actual_lengths() const
{
    return Real3(col_size() * HCP_X,
                 layer_size() * HCP_Y,
                 row_size() * voxel_radius() * 2);
}

// const Real volume() const
// {
//     return edge_lengths_[0] * edge_lengths_[1] * edge_lengths_[2];
// }

} // ecell4
