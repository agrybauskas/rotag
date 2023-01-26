#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/arginine-002.dump

$(dirname "$0")/../scripts/BondParameters_calculate_dihedral_angles ${pdbx_dump_file}
