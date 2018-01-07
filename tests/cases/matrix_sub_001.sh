#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

matrices_file=$(dirname "$0")/../inputs/matrices/matrices-008.dat

$(dirname "$0")/../scripts/matrix_sub ${matrices_file}
