%module LinearAlgebra
%{
  #include "LinearAlgebra.h"
%}

vector< vector<double> > create_ref_frame( double mid_atom_x,
					   double mid_atom_y,
					   double mid_atom_z,
					   double up_atom_x,
					   double up_atom_y,
					   double up_atom_z,
					   double side_atom_x,
					   double side_atom_y,
					   double side_atom_z );
