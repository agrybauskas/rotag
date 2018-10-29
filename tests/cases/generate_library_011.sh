#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/glutamine-H-rotation-only-002.dump
conf_model="rotation_only"
potential="composite"
energy_cutoff_atom=70.0
energy_cutoff_residue="Inf"
residue_id="2"
residue_chain="A"
pdbx_model="1"
residue_alt="."
small_angle="0.25*pi"
parameters='lj_epsilon=1.0,c_k=1.0,h_epsilon=1.0,r_sigma=2.0,cutoff_atom=0.5,cutoff_residue=1.0,cutoff_start=2.5,cutoff_end=5.0'

$(dirname "$0")/../scripts/generate_library ${residue_id} \
                                            ${residue_chain} \
                                            ${pdbx_model} \
                                            ${residue_alt} \
                                            ${conf_model} \
	                                    ${small_angle} \
	                                    ${potential} \
					    ${energy_cutoff_atom} \
					    ${energy_cutoff_residue} \
					    ${pdbx_dump_file} \
                                            ${parameters}
