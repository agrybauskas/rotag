#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/synthetic/xaa-rotation-only-002.dump
potential="hard_sphere"
interaction_specifier="label_atom_id N,C,CA,O"
energy_cutoff_atom=0
residue_id="18"
residue_chain="A"
pdbx_model_num="1"
residue_alt="."
angles="0.0..36.0..360.0"

$(dirname "$0")/../scripts/calc_favourable_angles \
               ${residue_id} \
               ${residue_chain} \
               ${pdbx_model_num} \
               ${residue_alt} \
	           ${angles} \
	           ${potential} \
               "${interaction_specifier}" \
	           ${energy_cutoff_atom} \
	           ${pdbx_dump_file}
