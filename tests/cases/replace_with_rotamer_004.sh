#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/lysine-rotation-only-001.dump
residue_unique_key="572,B,1,."
angle_values="chi1 0 & chi2 0.5*pi & chi3 pi & chi4 1.5*pi"

$(dirname "$0")/../scripts/replace_with_rotamer "${residue_unique_key}" \
	                                        "${angle_values}" \
					        ${pdbx_dump_file}
