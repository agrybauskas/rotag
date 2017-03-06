#!/bin/bash
cd "$(dirname "$0")"

symbols="x,y,z"
symb_matrices=../inputs/symb_matrices_004.dat

../programs/mult_matrix_product ${symbols} < ${symb_matrices}
