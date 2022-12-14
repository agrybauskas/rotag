#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/synthetic/single-atoms-001.dump
atom_i_id=1
atom_j_id=5

$(dirname "$0")/../scripts/coulomb ${atom_i_id} \
                                   ${atom_j_id} \
			           ${pdbx_dump_file}
