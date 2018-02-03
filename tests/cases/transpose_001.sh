#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

matrix_file=$(dirname "$0")/../inputs/matrices/matrix-001.dat

$(dirname "$0")/../scripts/transpose ${matrix_file}
