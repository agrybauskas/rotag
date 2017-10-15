#include <iostream>
#include <math.h>
#include <vector>

#include "LinearAlgebra.h"

using namespace std;

/* --------------------------------- Constants --------------------------------- */

/*
  Generates constants that are useful for other functions.
*/

/*
  Returns PI value.
  Input:
      none.
  Output:
      PI value.
*/

double pi()
{
    return 4 * atan2( 1, 1 );
}

/*
  Returns machine accuracy for floating point numbers.
  Input:
      none.
  Output:
      $EPSILON - machine accuracy value.
*/

double epsilon()
{
    double EPSILON = 1.0;

    while( ( 1.0 + 0.5 * EPSILON ) != 1.0 ) {
      EPSILON = 0.5 * EPSILON;
    }

    return EPSILON;
}

/* ------------------------- Numeric linear algebra -------------------------- */

/*
  Performs basic linear algebra operations for matrices.
 */

/*
  Creates local reference frame for any three given atoms positions in cartesian
  coordinate system.
  Input:
      {mid,up,side}_atom_{x,y,z} - Cartesian coordinates of three atoms.
  Output:
      local_ref_frame - Cartesian coordinates of points on x, y and z axis.
*/

vector< vector<double> > create_ref_frame( double mid_atom_x,
					   double mid_atom_y,
					   double mid_atom_z,
					   double up_atom_x,
					   double up_atom_y,
					   double up_atom_z,
					   double side_atom_x,
					   double side_atom_y,
					   double side_atom_z )
{
  /* Initializes 2D array. */
  vector< vector<double> > local_ref_frame;
  local_ref_frame.resize(3);
  for (int i = 0; i < 3; ++i) {
    local_ref_frame[i].resize(3);
  }

  /* Let local z-axis be colinear to bond between mid and up atoms. */
  local_ref_frame[2][0] = up_atom_x - mid_atom_x;
  local_ref_frame[2][1] = up_atom_y - mid_atom_y;
  local_ref_frame[2][2] = up_atom_z - mid_atom_z;

  /* Let local x-axis be perpendicular to mid-up and mid-side bonds. */
  local_ref_frame[0][0] =
    ( side_atom_y - mid_atom_y ) * local_ref_frame[2][2]
  - ( side_atom_z - mid_atom_z ) * local_ref_frame[2][1];
  local_ref_frame[0][1] =
  - ( side_atom_x - mid_atom_x ) * local_ref_frame[2][2]
  + ( side_atom_z - mid_atom_z ) * local_ref_frame[2][0];
  local_ref_frame[0][2] =
    ( side_atom_x - mid_atom_x ) * local_ref_frame[2][1]
  - ( side_atom_y - mid_atom_y ) * local_ref_frame[2][0];

  /* Let local y-axis be in the same plane as mid-up and mid-side bonds. */
  local_ref_frame[1][0] =
    local_ref_frame[2][1] * local_ref_frame[0][2]
  - local_ref_frame[2][2] * local_ref_frame[0][1];
  local_ref_frame[1][1] =
  - local_ref_frame[2][0] * local_ref_frame[0][2]
  + local_ref_frame[2][2] * local_ref_frame[0][0];
  local_ref_frame[1][2] =
    local_ref_frame[2][0] * local_ref_frame[0][1]
  - local_ref_frame[2][1] * local_ref_frame[0][0];

  return local_ref_frame;
}

/*
  Function calculates Euler rotational angles (alpha, beta, gamma) that are used
  to transform global reference frame to chosen one.
  Input  (1 arg): array of three atom coordinates in x, y, z form.
  Output (3 arg): euler angles (alpha, beta, gamma) in radians.
*/

vector<double> find_euler_angles( double mid_atom_x,
				  double mid_atom_y,
				  double mid_atom_z,
				  double up_atom_x,
				  double up_atom_y,
				  double up_atom_z,
				  double side_atom_x,
				  double side_atom_y,
				  double side_atom_z )
{
  vector<double> euler_angles;
  euler_angles.resize(3);

  double alpha_rad;
  double beta_rad;
  double gamma_rad;

  double z_axis_in_xy_plane;

  vector< vector<double> > local_ref_frame =
    create_ref_frame( mid_atom_x,  mid_atom_y,   mid_atom_z,
		      up_atom_x,   up_atom_y,    up_atom_z,
		      side_atom_x, side_atom_y,  side_atom_z );

  /* Projects local z-axis to global xy-plane. */
  z_axis_in_xy_plane =
    sqrt( local_ref_frame[2][0] * local_ref_frame[2][0]
	+ local_ref_frame[2][1] * local_ref_frame[2][1] );

  if( z_axis_in_xy_plane > epsilon() ) {
    alpha_rad =
      atan2( local_ref_frame[1][0] * local_ref_frame[2][1]
           - local_ref_frame[1][1] * local_ref_frame[2][0],
	     local_ref_frame[0][0] * local_ref_frame[2][1]
	   - local_ref_frame[0][1] * local_ref_frame[2][0] );
    beta_rad = atan2( z_axis_in_xy_plane, local_ref_frame[2][2] );
    gamma_rad =
      - atan2( - local_ref_frame[2][0], local_ref_frame[2][1] );
  } else {
    alpha_rad = 0.;
    beta_rad = ( local_ref_frame[2][2] > 0. ) ? 0. : pi();
    gamma_rad =
      - atan2( local_ref_frame[0][1], local_ref_frame[0][0] );
  }

  euler_angles[0] = alpha_rad;
  euler_angles[1] = beta_rad;
  euler_angles[2] = gamma_rad;
  
  return euler_angles;
}
