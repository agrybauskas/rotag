#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/synthetic/xaa-008.dump
residue_unique_key="18,A,1,."
angle_values="chi0 0 & chi1 pi"

$(dirname "$0")/../scripts/generate_rotamer "${residue_unique_key}" \
	                                    "${angle_values}" \
					    ${pdbx_dump_file}
