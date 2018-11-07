#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

unique_residue_key='315,A,1,.'
pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/arginine-002.dump

$(dirname "$0")/../scripts/identify_residue_atoms \
               ${unique_residue_key} \
               ${pdbx_dump_file}
