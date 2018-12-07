#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-rotation-only-001.dump
potential="hard_sphere"
interaction_specifier="label_atom_id N,C,CA,O"
energy_cutoff_atom=0
energy_cutoff_residue="Inf"
residue_id="18"
residue_chain="A"
pdbx_model_num="1"
residue_alt="."
small_angle="0.1*pi"

$(dirname "$0")/../scripts/calc_favourable_angles \
               ${residue_id} \
               ${residue_chain} \
               ${pdbx_model_num} \
               ${residue_alt} \
	           ${small_angle} \
	           ${potential} \
               "${interaction_specifier}" \
	           ${energy_cutoff_atom} \
               ${energy_cutoff_residue} \
	           ${pdbx_dump_file}
