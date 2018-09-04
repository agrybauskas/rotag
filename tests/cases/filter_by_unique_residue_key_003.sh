#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

unique_residue_key='.,O,7,.'
pdbx_dump_file=$(dirname "$0")/../inputs/5aox.dump

"$(dirname "$0")"/../scripts/filter_by_unique_residue_key ${unique_residue_key} \
                                                          ${pdbx_dump_file}
