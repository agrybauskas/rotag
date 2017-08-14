#!/bin/bash
cd "$(dirname "$0")"

vector_file=../inputs/two_3d_vectors_001.dat

../programs/vector_cross ${vector_file}
