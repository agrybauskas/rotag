#ifndef _LINEARALGEBRACPP_H_
#define _LINEARALGEBRACPP_H_

#include <vector>

std::vector< std::vector<double> > create_ref_frame( std::vector<double> mid_atom_coord,
                                                     std::vector<double> up_atom_coord,
                                                     std::vector<double> side_atom_coord );

double calc_vector_length( std::vector<double> vector );

#endif
