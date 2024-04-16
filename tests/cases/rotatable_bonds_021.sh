#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.dump
include_mainchain=0
include_hetatoms=1

"$(dirname "$0")"/../scripts/rotatable_bonds ${pdbx_dump_file} ${include_mainchain} ${include_hetatoms}
