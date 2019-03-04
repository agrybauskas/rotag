#include "LinearAlgebraCpp.h"
#include <iostream>
#include <math.h>
#include <vector>

/* ------------------------- Numeric linear algebra -------------------------- */

/*
  Creates local reference frame for any three given atoms positions in cartesian
  coordinate system.
  Input:
      {mid,up,side}_atom_coord - Cartesian coordinates of three atoms.
  Output:
      local_ref_frame - Cartesian coordinates of points on x, y and z axis.
*/

void create_ref_frame( std::vector<double> mid_atom_coord,
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

    // /* Normalizes all vectors to unit vectors. */
    // double vector_length = calculate_vector_length((double*) local_ref_frame[2]);
    // local_ref_frame[2][0] = local_ref_frame[2][0] / vector_length;
    // local_ref_frame[2][1] = local_ref_frame[2][1] / vector_length;
    // local_ref_frame[2][2] = local_ref_frame[2][2] / vector_length;

    // vector_length = calculate_vector_length((double*) local_ref_frame[0]);
    // local_ref_frame[0][0] = local_ref_frame[0][0] / vector_length;
    // local_ref_frame[0][1] = local_ref_frame[0][1] / vector_length;
    // local_ref_frame[0][2] = local_ref_frame[0][2] / vector_length;

    // vector_length = calculate_vector_length((double*) local_ref_frame[1]);
    // local_ref_frame[1][0] = local_ref_frame[1][0] / vector_length;
    // local_ref_frame[1][1] = local_ref_frame[1][1] / vector_length;
    // local_ref_frame[1][2] = local_ref_frame[1][2] / vector_length;

    // /* Converts to Perl data types. */
    // AV* local_ref_frame_av = new AV();
    // for( int i = 0; i < 3; i++ ) {
    //     for( int j = 0; j < 3; j++ ) {

    //     }
    // }
}

// /*
//   Calculates vector length.
//   Input:
//       vector - 1x3 (if 1x4, last column is ignored) matrix.
//   Output:
//       vector length.
// */

// double calculate_vector_length( double* vector_ptr )
// {
//     double vector[3] = {vector_ptr[0], vector_ptr[1], vector_ptr[2]};
//     return sqrt(pow( vector[0], 2 ) + pow( vector[1], 2 ) + pow( vector[2], 2 ));
// }
