#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

atom_coords_file=$(dirname "$0")/../inputs/atom-coords-004.dat

$(dirname "$0")/../scripts/bond_angle ${atom_coords_file}
