#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-rotation-only-001.dump
conf_model="rotation_only"
potential="composite"
energy_cutoff_atom=1.0
residue_id="18"
residue_chain="A"
residue_entity="1"
residue_alt="."
angles="18.0"
parameters='lj_k=1.0,c_k=1.0,h_k=1.0,t_k=1.0,cutoff_atom=0.5,cutoff_start=2.5,cutoff_end=5.0'

$(dirname "$0")/../scripts/generate_library ${residue_id} \
                                            ${residue_chain} \
                                            ${residue_entity} \
                                            ${residue_alt} \
                                            ${conf_model} \
	                                        ${angles} \
	                                        ${potential} \
					                        ${energy_cutoff_atom} \
					                        ${pdbx_dump_file} \
                                            ${parameters}
