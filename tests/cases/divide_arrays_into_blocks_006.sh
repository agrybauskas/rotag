#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

vector_file=$(dirname "$0")/../inputs/matrices/vector-003.dat
threads=3

$(dirname "$0")/../scripts/divide_arrays_into_blocks ${vector_file} ${threads}
