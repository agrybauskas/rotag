#!/bin/bash
cd "$(dirname "$0")"

symbols=""
symb_matrices=../inputs/symb_matrices_002.dat

../programs/two_matrix_product ${symbols} < ${symb_matrices}
