#!/bin/bash
cd "$(dirname "$0")"

symbols=""
symb_matrices=../input/symb_matrices_002.dat

../programs/symb_dot_product ${symbols} < ${symb_matrices}
