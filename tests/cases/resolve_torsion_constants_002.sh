#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.dump
atom_ids="1414,1413,1412,1409"

$(dirname "$0")/../scripts/resolve_torsion_constants ${atom_ids} ${pdbx_dump_file} 2>&1
