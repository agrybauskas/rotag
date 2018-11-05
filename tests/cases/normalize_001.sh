#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

vector_file=$(dirname "$0")/../inputs/matrices/vector-001.dat

$(dirname "$0")/../scripts/normalize ${vector_file}
