#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/threonine-rotation-only-001.dump
residue_id="340"
angle_values="chi0 0"

$(dirname "$0")/../scripts/generate_rotamer "${residue_id}" \
	                                    "${angle_values}" \
					    ${pdbx_dump_file}
