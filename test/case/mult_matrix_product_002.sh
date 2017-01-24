#!/bin/bash
cd "$(dirname "$0")"

symbols="x,y,z,chi"
symb_matrices=../input/symb_matrices_005.dat

../programs/mult_matrix_product ${symbols} < ${symb_matrices}
