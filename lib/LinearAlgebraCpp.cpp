#include "LinearAlgebraCpp.h"

#include <map>
#include <iostream>
#include <limits>
#include <math.h>
#include <vector>

#include "AlgebraicMatrix.h"

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
  double vector_length;
  for ( int i = 0; i < 3; i++ ) {
    vector_length = calc_vector_length( local_ref_frame[i] );
    for ( int j = 0; j < 3; j++ ) {
      local_ref_frame[i][j] = local_ref_frame[i][j] / vector_length;
    }
  }

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

  if ( z_axis_in_xy_plane > std::numeric_limits<double>::epsilon() ) {
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

AlgebraicMatrix matrix_product( AlgebraicMatrix left_matrix,
                                AlgebraicMatrix right_matrix,
                                std::map<std::string, double> symbol_values )
{
  /* First, evaluates all matrices if they are not evaluated. */
  if ( left_matrix.get_is_evaluated() != 1 ) {
    left_matrix.evaluate( symbol_values );
  }
  if ( right_matrix.get_is_evaluated() != 1 ) {
    right_matrix.evaluate( symbol_values );
  }

  std::vector< std::vector<double> > local_left_matrix = left_matrix.get_matrix();
  std::vector< std::vector<double> > local_right_matrix = right_matrix.get_matrix();

  /* Notifies error, when the column number of left matrix does not equal the
     row number of the right matrix. */
  if ( local_left_matrix[0].size() != local_right_matrix.size() ) {
    std::cout << "A row number of a left matrix is NOT equal to the column\n"
              << "number of the right matrix."
              << std::endl;
    exit( EXIT_FAILURE );
  }

  std::vector< std::vector<double> > local_matrix_product(
      local_left_matrix.size(), std::vector<double>( local_right_matrix[0].size() ) );

  for ( int left_row = 0; left_row < local_left_matrix.size(); left_row++  ) {
    for ( int right_col = 0; right_col < local_right_matrix[0].size(); right_col++ ) {
      for ( int right_row = 0; right_row < local_right_matrix.size(); right_row++ ) {
        local_matrix_product[left_row][right_col] +=
          local_left_matrix[left_row][right_row] *
          local_right_matrix[right_row][right_col];
      }
    }
  }

  AlgebraicMatrix matrix_product_obj;
  matrix_product_obj.set_matrix( local_matrix_product );

  return matrix_product_obj;
}

void mult_matrix_product( std::vector<AlgebraicMatrix> matrices,
                          std::map<std::string, double> symbol_values )
{
  std::vector<AlgebraicMatrix> mult_matrix_product;

  if ( matrices.size() == 1 ) {
    AlgebraicMatrix matrix = matrices[0];
    if ( ! matrix.get_is_evaluated() ) {
      matrix.evaluate( symbol_values );
    }
    mult_matrix_product.push_back(matrix);
  } else {

  }
}
