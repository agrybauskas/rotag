#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/lysine-rotation-only-001.dump
conf_model="rotation_only"
potential="composite"
residue_id="572"
residue_chain="B"
residue_entity="1"
residue_alt="."
angles="..180.0.."
parameters='lj_k=1.0,c_k=1.0,h_k=1.0,t_k=1.0,cutoff_atom=70.0,cutoff_start=2.5,cutoff_end=5.0'
threads="1"
auth_residue_id="617"
auth_residue_chain="B"

$(dirname "$0")/../scripts/generate_library ${residue_id} \
                                            ${residue_chain} \
                                            ${residue_entity} \
                                            ${residue_alt} \
                                            ${conf_model} \
	                                        ${angles} \
	                                        ${potential} \
					                        ${pdbx_dump_file} \
                                            ${parameters} \
                                            ${threads} \
                                            ${auth_residue_id} \
                                            ${auth_residue_chain}
