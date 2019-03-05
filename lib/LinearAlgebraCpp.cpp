#include "LinearAlgebraCpp.h"

#include <iostream>
#include <math.h>
#include <vector>

/* ------------------------- Numeric linear algebra -------------------------- */

std::vector< std::vector<double> > create_ref_frame( std::vector<double> mid_atom_coord,
                                                     std::vector<double> up_atom_coord,
                                                     std::vector<double> side_atom_coord )
{
    std::vector< std::vector<double> > local_ref_frame(3, std::vector<double>(4) );

    /* Let local z-axis be colinear to bond between mid and up atoms. */
    local_ref_frame[2][0] = up_atom_coord[0] - mid_atom_coord[0];
    local_ref_frame[2][1] = up_atom_coord[1] - mid_atom_coord[1];
    local_ref_frame[2][2] = up_atom_coord[2] - mid_atom_coord[2];

    /* Let local x-axis be perpendicular to mid-up and mid-side bonds. */
    local_ref_frame[0][0] =
        ( side_atom_coord[1] - mid_atom_coord[1] ) * local_ref_frame[2][2] -
        ( side_atom_coord[2] - mid_atom_coord[2] ) * local_ref_frame[2][1];
    local_ref_frame[0][1] =
      - ( side_atom_coord[0] - mid_atom_coord[0] ) * local_ref_frame[2][2] +
        ( side_atom_coord[2] - mid_atom_coord[2] ) * local_ref_frame[2][0];
    local_ref_frame[0][2] =
        ( side_atom_coord[0] - mid_atom_coord[0] ) * local_ref_frame[2][1] -
        ( side_atom_coord[1] - mid_atom_coord[1] ) * local_ref_frame[2][0];

    /* Let local y-axis be in the same plane as mid-up and mid-side bonds. */
    local_ref_frame[1][0] =
        local_ref_frame[2][1] * local_ref_frame[0][2] -
        local_ref_frame[2][2] * local_ref_frame[0][1];
    local_ref_frame[1][1] =
      - local_ref_frame[2][0] * local_ref_frame[0][2] +
        local_ref_frame[2][2] * local_ref_frame[0][0];
    local_ref_frame[1][2] =
        local_ref_frame[2][0] * local_ref_frame[0][1] -
        local_ref_frame[2][1] * local_ref_frame[0][0];

    /* Normalizes all vectors to unit vectors. */
    double vector_length = calc_vector_length( local_ref_frame[2] );
    local_ref_frame[2][0] = local_ref_frame[2][0] / vector_length;
    local_ref_frame[2][1] = local_ref_frame[2][1] / vector_length;
    local_ref_frame[2][2] = local_ref_frame[2][2] / vector_length;

    vector_length = calc_vector_length( local_ref_frame[0] );
    local_ref_frame[0][0] = local_ref_frame[0][0] / vector_length;
    local_ref_frame[0][1] = local_ref_frame[0][1] / vector_length;
    local_ref_frame[0][2] = local_ref_frame[0][2] / vector_length;

    vector_length = calc_vector_length( local_ref_frame[1] );
    local_ref_frame[1][0] = local_ref_frame[1][0] / vector_length;
    local_ref_frame[1][1] = local_ref_frame[1][1] / vector_length;
    local_ref_frame[1][2] = local_ref_frame[1][2] / vector_length;

    return local_ref_frame;
}

double calc_vector_length( std::vector<double> vector )
{
    return sqrt( pow( vector[0], 2 ) +
                 pow( vector[1], 2 ) +
                 pow( vector[2], 2 ) );
}
