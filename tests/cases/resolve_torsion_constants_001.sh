#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-001.dump
atom_ids="152,151,148,147"

$(dirname "$0")/../scripts/resolve_torsion_constants ${atom_ids} ${pdbx_dump_file} 2>&1
