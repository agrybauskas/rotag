#!/bin/bash
cd "$(dirname "$0")"

symbols="x,y,z"
symb_matrices=../input/symb_matrices_001.dat

../programs/two_matrix_product ${symbols} < ${symb_matrices}
