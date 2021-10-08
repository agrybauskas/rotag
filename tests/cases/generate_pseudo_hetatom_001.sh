#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-rotation-translation-001.dump
atom_id=152
angle_ranges="chi1 -1*pi,1*pi & r1 -0.1,0.1 & r2 -0.1,0.1 & theta1 -0.1*pi.0.1*pi & theta2 -0.1*pi,0.1*pi"
num_of_angles="chi1 20 & r1 1 & r2 1 & theta1 1 & theta2 1"

$(dirname "$0")/../scripts/generate_pseudo_hetatom "${atom_id}" \
                                                   "${angle_ranges}" \
                                                   "${num_of_angles}" \
                                                   ${pdbx_dump_file}
