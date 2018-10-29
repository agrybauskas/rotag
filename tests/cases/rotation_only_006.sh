#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/synthetic/xaa-002.dump
atom_id=18

$(dirname "$0")/../scripts/rotation_only ${atom_id} ${pdbx_dump_file}
