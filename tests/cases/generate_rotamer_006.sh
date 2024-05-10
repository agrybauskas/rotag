#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/glutamic-acid-rotation-only-001.dump
residue_unique_key="14,A,1,.,59,A"
angle_values="chi1 0 & chi2 pi & chi3 0"

$(dirname "$0")/../scripts/generate_rotamer "${residue_unique_key}" \
	                                    "${angle_values}" \
					    ${pdbx_dump_file}
