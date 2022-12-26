#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-with-connections-006.dump
atom_id=1925
angle_and_length_ranges="N-CA 0.5*pi,1.0*pi & CA-MG 0.0,2.0"
num_of_angles_and_lengths="N-CA 4 & CA-MG 3"
no_full_range=1

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file} \
                                                   ${no_full_range}
