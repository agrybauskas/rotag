#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/histidine-001.dump

"$(dirname "$0")"/../scripts/stretchable_bonds ${pdbx_dump_file}
