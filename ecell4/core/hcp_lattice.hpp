#ifndef __ECELL4_HCP_LATTICE_HPP
#define __ECELL4_HCP_LATTICE_HPP

#include "types.hpp"
#include "Integer3.hpp"
#include "Real3.hpp"

namespace ecell4
{

inline Real3
calc_lattice_position(const Real& voxel_radius, const Integer3& point)
throw ()
{
    const Integer &col(point[0]),
                  &row(point[1]),
                  &layer(point[2]);
    return Real3(col * sqrt(8.0 / 3.0),
                 (col % 2) / sqrt(3.0) + layer * sqrt(3.0),
                 row * 2 + (layer + col) % 2)
           * voxel_radius;
}

inline Integer
calc_neighbor(const Integer3& size, const Integer& coord, const Integer& idx)
throw (NotFound)
{
    const Integer &num_cols(size[0]),
                  &num_rows(size[1]);
    const Integer layer_size(num_cols * num_rows);
    const Integer odd_col(((coord % layer_size) / num_rows) & 1),
                  odd_layer((coord / layer_size) & 1);

    switch (idx)
    {
        case 0:
            return coord - 1;
        case 1:
            return coord + 1;
        case 2:
            return coord + (odd_col ^ odd_layer) - num_rows - 1;
        case 3:
            return coord + (odd_col ^ odd_layer) - num_rows;
        case 4:
            return coord + (odd_col ^ odd_layer) + num_rows - 1;
        case 5:
            return coord + (odd_col ^ odd_layer) + num_rows;
        case 6:
            return coord - (2 * odd_col - 1) * layer_size - num_rows;
        case 7:
            return coord - (2 * odd_col - 1) * layer_size + num_rows;
        case 8:
            return coord + (odd_col ^ odd_layer) - layer_size - 1;
        case 9:
            return coord + (odd_col ^ odd_layer) - layer_size;
        case 10:
            return coord + (odd_col ^ odd_layer) + layer_size - 1;
        case 11:
            return coord + (odd_col ^ odd_layer) + layer_size;
    }

    throw NotFound("Invalid argument: idx");
}

} // ecell4

#endif /* __ECELL4_HCP_LATTICE_HPP */
