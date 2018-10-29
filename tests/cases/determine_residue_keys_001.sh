#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-002.dump

$(dirname "$0")/../scripts/determine_residue_keys 1 ${pdbx_file}
