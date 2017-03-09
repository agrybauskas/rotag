#!/bin/bash
cd "$(dirname "$0")"

two_vectors=../inputs/two_3d_vectors_001.dat

../programs/vector_cross < ${two_vectors}
