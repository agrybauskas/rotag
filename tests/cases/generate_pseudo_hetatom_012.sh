#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-008.dump
atom_id=152
angle_and_length_ranges="CA-CB-OG 0.5*pi,1.0*pi & CA-CB 1.0,1.5 & chi1 0.0,0.5*pi"
num_of_angles_and_lengths="CA-CB-OG 3 & CA-CB 3 & chi1 3"
full_range="CA-CB-OG 1 & CA-CB 1 & chi1 0"

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file} \
                                                   "${full_range}"
