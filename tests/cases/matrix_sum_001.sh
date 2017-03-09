#!/bin/bash
cd "$(dirname "$0")"

two_vectors=../inputs/two_matrices_001.dat

../programs/matrix_sum < ${two_vectors}
