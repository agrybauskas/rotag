#!/bin/bash
cd "$(dirname "$0")"

symbols="chi,psi,theta"
symb_matrices=../input/symb_matrices_005.dat

../programs/mult_matrix_product ${symbols} < ${symb_matrices}
