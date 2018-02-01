#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/cysteine-001.dump

"$(dirname "$0")"/../scripts/rotatable_bonds ${pdbx_dump_file}
