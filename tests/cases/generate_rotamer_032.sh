#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/lysine-H-rotation-only-001.dump
residue_unique_key="572,B,1,."
angle_values="chi0 0 & chi1 pi & chi2 0 & chi3 pi & chi4 pi"

$(dirname "$0")/../scripts/generate_rotamer "${residue_unique_key}" \
	                                    "${angle_values}" \
					    ${pdbx_dump_file}
