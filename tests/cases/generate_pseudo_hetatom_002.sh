#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-bond-stretching-only-001.dump
atom_id=152
angle_and_length_ranges="r1 -0.1,0.1 & r2 -0.1,0.1"
num_of_angles_and_lengths="r1 2 & r2 2"

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file}
