#!/bin/bash
cd "$(dirname "$0")"

symb_matrices=../input/symb_matrices_001.dat

../programs/two_matrix_product < ${symb_matrices}
