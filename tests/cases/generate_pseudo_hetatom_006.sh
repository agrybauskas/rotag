#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-with-connections-008.dump
atom_id=1925
angle_and_length_ranges="N-CA-MG 0.5*pi,1.0*pi & CA-MG 0.5,1.0"
num_of_angles_and_lengths="N-CA-MG 5 & CA-MG 5"

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file}
