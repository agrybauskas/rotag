#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

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

double create_ref_frame( double mid_atom_x,
			 double mid_atom_y,
			 double mid_atom_z,
			 double up_atom_x,
			 double up_atom_y,
			 double up_atom_z,
			 double side_atom_x,
			 double side_atom_y,
			 double side_atom_z )
{
  double local_ref_frame[3][3];

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
}

main()
{
  return 0;
}
