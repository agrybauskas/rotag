#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/lysine-rotation-only-001.dump
residue_unique_key="572,B,1,.,617,B"
angle_values="chi1 0 & chi2 pi & chi3 0 & chi4 pi"

$(dirname "$0")/../scripts/generate_rotamer "${residue_unique_key}" \
	                                    "${angle_values}" \
					    ${pdbx_dump_file}
