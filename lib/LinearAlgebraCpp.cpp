#include "LinearAlgebraCpp.h"
#include <iostream>
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

/* ------------------------- Numeric linear algebra -------------------------- */

/*
  Creates local reference frame for any three given atoms positions in cartesian
  coordinate system.
  Input:
      {mid,up,side}_atom_coord - Cartesian coordinates of three atoms.
  Output:
      local_ref_frame - Cartesian coordinates of points on x, y and z axis.
*/

void create_ref_frame()
{
    std::cout << "Haha\n" << std::endl;
}
