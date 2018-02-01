#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/lysine-001.dump
atom_id=8649

$(dirname "$0")/../scripts/rotation_only ${atom_id} ${pdbx_dump_file}
