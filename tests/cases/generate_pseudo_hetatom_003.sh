#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-angle-bending-only-001.dump
atom_id=152
angle_and_length_ranges="theta1 -0.1*pi,0.1*pi & theta2 -0.1*pi,0.1*pi"
num_of_angles_and_lengths="theta1 2 & theta2 2"

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file}
