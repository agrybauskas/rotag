#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-007.dump
atom_id=152
angle_and_length_ranges="N-CA-CB 0.5*pi,1.0*pi & CA-CB-OG 0.5*pi,1.0*pi"
num_of_angles_and_lengths="N-CA-CB 5 & CA-CB-OG 5"
no_full_range=1

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file} \
                                                   ${no_full_range}
