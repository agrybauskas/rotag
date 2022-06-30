#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-with-connections-002.dump
atom_id=1925

$(dirname "$0")/../scripts/rotation_translation ${atom_id} ${pdbx_dump_file}
