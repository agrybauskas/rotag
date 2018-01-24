#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-rotation-only-001.dump
atom_id=1414
angle_ranges="chi0 0,2*pi & chi1 0,2*pi"
num_of_angles="chi0 4 & chi1 4"

$(dirname "$0")/../scripts/generate_pseudo "${atom_id}" \
	                                   "${angle_ranges}" \
	                                   "${num_of_angles}" \
					   ${pdbx_dump_file}
