#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

atom_coords_file=$(dirname "$0")/../inputs/atom-coords-002.dat

$(dirname "$0")/../scripts/rmsd ${atom_coords_file}
