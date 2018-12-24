#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/lysine-rotation-only-001.dump
atom_id=8649
angle_ranges="chi1 0,2*pi & chi2 0,2*pi & chi3 0,2*pi & chi4 0,2*pi"
num_of_angles="chi1 3 & chi2 3 & chi3 3 & chi4 3"

$(dirname "$0")/../scripts/generate_pseudo "${atom_id}" \
                                           "${angle_ranges}" \
                                           "${num_of_angles}" \
                                           ${pdbx_dump_file}
