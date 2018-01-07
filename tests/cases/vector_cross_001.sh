#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

vectors_3d_file=$(dirname "$0")/../inputs/matrices/vectors-3d-001.dat

$(dirname "$0")/../scripts/vector_cross ${vectors_3d_file}
