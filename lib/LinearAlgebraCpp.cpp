#include "LinearAlgebraCpp.h"

#include <iostream>
#include <limits>
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

std::vector<double> find_euler_angle( std::vector<double> mid_atom_coord,
                                      std::vector<double> up_atom_coord,
                                      std::vector<double> side_atom_coord )
{
  double alpha_rad;
  double beta_rad;
  double gamma_rad;

  double z_axis_in_xy_plane;

  std::vector< std::vector<double> > local_ref_frame =
    create_ref_frame( mid_atom_coord, up_atom_coord, side_atom_coord );

  /* Projects local z-axis to global xy-plane. */
  z_axis_in_xy_plane =
     sqrt( local_ref_frame[2][0] * local_ref_frame[2][0] +
           local_ref_frame[2][1] * local_ref_frame[2][1] );

  if( z_axis_in_xy_plane > std::numeric_limits<double>::epsilon() ) {
    alpha_rad =
      atan2( local_ref_frame[1][0] * local_ref_frame[2][1] -
             local_ref_frame[1][1] * local_ref_frame[2][0],
             local_ref_frame[0][0] * local_ref_frame[2][1] -
             local_ref_frame[0][1] * local_ref_frame[2][0] );
    beta_rad = atan2( z_axis_in_xy_plane, local_ref_frame[2][2] );
    gamma_rad = - atan2( - local_ref_frame[2][0], local_ref_frame[2][1] );
  } else {
    alpha_rad = 0.;
    beta_rad = ( local_ref_frame[2][2] > 0. ) ? 0. : M_PI;
    gamma_rad = - atan2( local_ref_frame[0][1], local_ref_frame[0][0] );
  }

  return std::vector<double> { alpha_rad, beta_rad, gamma_rad };
}

double calc_vector_length( std::vector<double> vector )
{
  return sqrt( pow( vector[0], 2 ) + pow( vector[1], 2 ) + pow( vector[2], 2 ) );
}

/* ------------------------ Symbolic linear algebra -------------------------- */

std::vector< std::vector<double> > transpose( std::vector< std::vector<double> > matrix )
{
  int row_count = matrix.size();
  int col_count = matrix[0].size();

  std::vector< std::vector<double> > transposed_matrix ( row_count,
                                                         std::vector<double> ( col_count ) );

  for ( int row = 0; row < transposed_matrix.size(); row++ ) {
    for ( int col = 0; col < transposed_matrix[row].size(); col++ ) {
      transposed_matrix[col][row] = matrix[row][col];
    }
  }

  return transposed_matrix;
}
