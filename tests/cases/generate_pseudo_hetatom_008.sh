#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-004.dump
atom_id=152
angle_and_length_ranges="r1 1.0,2.0 & r2 1.0,2.0"
num_of_angles_and_lengths="r1 5 & r2 5"
no_full_range=1

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_and_length_ranges}" \
                                                   "${num_of_angles_and_lengths}" \
                                                   ${pdbx_dump_file} \
                                                   ${no_full_range}
