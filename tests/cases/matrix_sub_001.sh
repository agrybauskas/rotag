#!/bin/bash
cd "$(dirname "$0")"

matrix_file=../inputs/two_matrices_002.dat

../programs/matrix_sub ${matrix_file}
