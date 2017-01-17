#!/bin/bash
cd "$(dirname "$0")"

symb_matrices=../input/symb_matrices_005.dat

../programs/mult_matrix_product < ${symb_matrices}
