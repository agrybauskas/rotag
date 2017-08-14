#!/bin/bash
cd "$(dirname "$0")"

matrix_file=../inputs/two_matrices_001.dat

../programs/matrix_sum ${matrix_file}
