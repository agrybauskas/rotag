#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-with-connections-009.dump
atom_id=1925
angle_and_length_ranges="tau1 0.5*pi,1.0*pi & r1 0.5,1.0 & theta1 0.5*pi,1.0*pi"
num_of_angles_and_lengths="tau1 5 & r1 5 & theta1 5"
no_full_range=1

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file} \
                                                   ${no_full_range}
