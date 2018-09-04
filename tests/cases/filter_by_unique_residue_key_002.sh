#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

unique_residue_key='81,B,1,.'
pdbx_dump_file=$(dirname "$0")/../inputs/1bvi.dump

"$(dirname "$0")"/../scripts/filter_by_unique_residue_key ${unique_residue_key} \
                                                          ${pdbx_dump_file}
