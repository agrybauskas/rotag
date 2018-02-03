#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

set_file=$(dirname "$0")/../inputs/sets-of-angles-001.dat

$(dirname "$0")/../scripts/permutation ${set_file}
