#!/bin/bash
cd "$(dirname "$0")"

symb_matrices=../input/symb_matrices_001.dat

../programs/symb_dot_product < ${symb_matrices}
