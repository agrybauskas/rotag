#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-rotation-translation-001.dump
atom_id=152
angle_and_length_ranges="chi1 -pi,pi & r1 1,2 & r2 1,2 & theta1 -0.1*pi,0.1*pi & theta2 -pi,pi"
num_of_angles_and_lengths="chi1 20 & r1 5 & r2 5 & theta1 20 & theta2 20"

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file}
