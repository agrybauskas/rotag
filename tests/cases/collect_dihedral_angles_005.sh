#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-H-001.dump

$(dirname "$0")/../scripts/collect_dihedral_angles ${pdbx_dump_file}
