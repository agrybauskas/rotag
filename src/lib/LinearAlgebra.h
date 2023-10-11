#ifndef SRC_LIB_LINEARALGEBRA_H_
#define SRC_LIB_LINEARALGEBRA_H_

#include <map>
#include <string>
#include <vector>

#include "AlgebraicMatrix.h"

// ------------------------- Numeric linear algebra ----------------------------

/*
  Creates local reference frame for any three given atoms positions in cartesian
  coordinate system.
  Input:
      {mid,up,side}_atom_coord - Cartesian coordinates of three atoms.
  Output:
      local_ref_frame - Cartesian coordinates of points on x, y and z axis.
*/

std::vector<std::vector<double>>
create_ref_frame(std::vector<double> mid_atom_coord,
                 std::vector<double> up_atom_coord,
                 std::vector<double> side_atom_coord);

/*
  Function calculates Euler rotational angles (alpha, beta, gamma) that are used
  to transform global reference frame to chosen one.
  Input:
      {mid,up,side}_atom_coord - Cartesian coordinates of three atoms.
  Output:
      euler angles (alpha, beta, gamma) in radians.
*/

std::vector<double> find_euler_angle(std::vector<double> mid_atom_coord,
                                     std::vector<double> up_atom_coord,
                                     std::vector<double> side_atom_coord);

/*
  Calculates vector length.
  Input:
      vector - 1x3 (if 1x4, last column is ignored) matrix.
  Output:
      vector length.
*/

double calc_vector_length(std::vector<double> vector);

// ------------------------ Symbolic linear algebra ----------------------------

/*
  Performs basic linear algebra on symbolic expressions.

  Example of rotation along z-axis by chi angle:

       / cos(chi) -sin(chi) 0 \   / x \   / x * cos(chi) + y * sin(chi) \
       | sin(chi)  cos(chi) 0 | * | y | = | x * sin(chi) + y * cos(chi) |
       \    0         0     1 /   \ z /   \              z              /
*/

/*
  Transposes matrix.
  Input:
      matrix - array representing matrix.
  Output:
      transposed_matrix - transposed matrix.
*/

std::vector<std::vector<double>>
transpose(std::vector<std::vector<double>> matrix);

/*
  Calculates matrix product of two matrices.
  Input:
      {left, right}_matrix - matrices.
      symbol_values - values of the unknown variable(-s).
  Output:
      matrix_product - matrix product.
*/

AlgebraicMatrix matrix_product(AlgebraicMatrix left_matrix,
                               AlgebraicMatrix right_matrix,
                               std::map<std::string, double> symbol_values);

/*
  Calculates matrix product of list of any size of matrices.
  Input:
      matrices - list of matrices.
      symbol_values - values of the unknown variable(-s).
  Output:
      mult_matrix_product - matrix product.
*/

std::vector<AlgebraicMatrix>
mult_matrix_product(std::vector<AlgebraicMatrix> matrices,
                    std::map<std::string, double> symbol_values);

#endif  // SRC_LIB_LINEARALGEBRA_H_
