#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-rotation-only-001.dump
residue_unique_key="18,A,1,.,18,A"
angle_values="chi1 0.5*pi"

$(dirname "$0")/../scripts/replace_with_rotamer \
               "${residue_unique_key}" \
     	       "${angle_values}" \
		       ${pdbx_dump_file}
