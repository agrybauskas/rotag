#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

atom_coords_file=$(dirname "$0")/../inputs/atom-coords-003.dat

$(dirname "$0")/../scripts/bond_length ${atom_coords_file}
