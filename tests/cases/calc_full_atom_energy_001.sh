#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-rotation-only-001.dump
potential="hard_sphere"
interaction_specifier="label_atom_id N,C,CA,O"
energy_cutoff_atom=0
residue_id="18"
residue_chain="A"
pdbx_model_num="1"
residue_alt="."
checkable_angles="1.0*pi;1.1*pi;1.2*pi;1.3*pi;1.4*pi;1.5*pi;1.6*pi;1.7*pi;1.8*pi;1.9*pi;2.0*pi"

$(dirname "$0")/../scripts/calc_full_atom_energy \
               ${residue_id} \
               ${residue_chain} \
               ${pdbx_model_num} \
               ${residue_alt} \
	           ${potential} \
               "${interaction_specifier}" \
	           ${energy_cutoff_atom} \
               ${checkable_angles} \
	           ${pdbx_dump_file}
