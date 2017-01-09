#!/bin/bash
cd "$(dirname "$0")"

symbols="x,y,z"
rec_symb_matrices=../input/rec_symb_matrices_001.dat

../programs/rec_symb_dot_product ${symbols} < ${rec_symb_matrices}
