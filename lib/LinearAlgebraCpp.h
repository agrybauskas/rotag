#ifndef _LINEARALGEBRACPP_H_
#define _LINEARALGEBRACPP_H_

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

std::vector< std::vector<double> > create_ref_frame( std::vector<double> mid_atom_coord,
                                                     std::vector<double> up_atom_coord,
                                                     std::vector<double> side_atom_coord );

/*
  Function calculates Euler rotational angles (alpha, beta, gamma) that are used
  to transform global reference frame to chosen one.
  Input:
      {mid,up,side}_atom_coord - Cartesian coordinates of three atoms.
  Output:
      euler angles (alpha, beta, gamma) in radians.
*/

std::vector<double> find_euler_angle( std::vector<double> mid_atom_coord,
                                      std::vector<double> up_atom_coord,
                                      std::vector<double> side_atom_coord );

/*
  Calculates vector length.
  Input:
      vector - 1x3 (if 1x4, last column is ignored) matrix.
  Output:
      vector length.
*/

double calc_vector_length( std::vector<double> vector );

#endif
