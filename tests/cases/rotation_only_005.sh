#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/arginine-002.dump
atom_id=2179

$(dirname "$0")/../scripts/rotation_only ${atom_id} ${pdbx_dump_file}
