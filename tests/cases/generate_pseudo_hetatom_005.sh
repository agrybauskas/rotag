#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-with-connections-007.dump
atom_id=1925
angle_and_length_ranges="CB-CG-OD1.MG 0.5*pi,1.0*pi & CG-OD1.MG 0.5*pi,1.0*pi"
num_of_angles_and_lengths="CB-CG-OD1.MG 5 & CG-OD1.MG 5"
full_range="CB-CG-OD1.MG 0 & CG-OD1.MG 1"

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file} \
                                                   "${full_range}"
