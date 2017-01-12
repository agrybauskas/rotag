#!/bin/bash
cd "$(dirname "$0")"

symbols=""
symb_matrices=../input/symb_matrices_003.dat

../programs/two_matrix_product ${symbols} < ${symb_matrices}
