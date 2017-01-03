#!/bin/bash
cd "$(dirname "$0")"

symbols="x,y,z"
symb_matrices=../input/symb_matrices_001.dat

../programs/symb_dot_product ${symbols} < ${symb_matrices}
