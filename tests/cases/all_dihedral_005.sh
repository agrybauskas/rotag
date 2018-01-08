#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-H-001.dump

$(dirname "$0")/../scripts/all_dihedral ${pdbx_dump_file}
