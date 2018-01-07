#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

matrices_file=$(dirname "$0")/../inputs/matrices/matrices-007.dat

$(dirname "$0")/../scripts/matrix_sum ${matrices_file}
