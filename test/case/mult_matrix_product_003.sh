#!/bin/bash
cd "$(dirname "$0")"

symbols="x,y,z,chi,i"
symb_matrices=../input/symb_matrices_006.dat

../programs/mult_matrix_product ${symbols} < ${symb_matrices}
