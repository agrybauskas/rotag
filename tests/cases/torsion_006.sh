#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.dump
atom_i_id=1414
debug=1

$(dirname "$0")/../scripts/torsion ${atom_i_id} ${pdbx_dump_file} ${debug} 2>&1
