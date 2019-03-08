#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-001.dump
atom_i_id=152

$(dirname "$0")/../scripts/torsion_new ${atom_i_id} ${pdbx_dump_file}
