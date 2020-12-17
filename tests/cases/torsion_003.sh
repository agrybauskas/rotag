#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-003.dump
atom_i_id=152

$(dirname "$0")/../scripts/torsion ${atom_i_id} ${pdbx_dump_file}
