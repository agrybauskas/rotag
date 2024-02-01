#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-with-connections-009.dump
atom_id=1925
angle_and_length_ranges="CB-CG-OD1.MG 0.5*pi,1.0*pi & OD1.MG 0.5,2.5 & CG-OD1.MG 0.5*pi,1.0*pi"
num_of_angles_and_lengths="CB-CG-OD1.MG 4 & OD1.MG 4 & CG-OD1.MG 4"
full_range="CB-CG-OD1.MG 0 & OD1.MG 1 & CG-OD1.MG 1"

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file} \
                                                   "${no_full_range}"
