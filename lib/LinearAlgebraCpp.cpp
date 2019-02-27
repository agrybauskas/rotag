#include "LinearAlgebraCpp.h"
#include <iostream>
#include <string>
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

void create_ref_frame( SV *mid_atom_coord_ptr, SV *up_atom_coord_ptr, SV *side_atom_coord_ptr )
{
    AV *mid_atom_coord_av = (AV*) SvRV( mid_atom_coord_ptr );
    AV *up_atom_coord_av = (AV*) SvRV( up_atom_coord_ptr );
    AV *side_atom_coord_av = (AV*) SvRV( side_atom_coord_ptr );

    double mid_atom_coord[3];
    double up_atom_coord[3];
    double side_atom_coord[3];
}
