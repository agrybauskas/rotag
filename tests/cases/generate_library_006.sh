#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/glutamic-acid-001.dump
conf_model="rotation_only"
potential="composite"
energy_cutoff_atom=1.0
energy_cutoff_residue="Inf"
residue_id="14"
residue_chain="A"
residue_entity="1"
residue_alt="."
small_angle="0.5*pi"

$(dirname "$0")/../scripts/generate_library ${residue_id} \
                                            ${residue_chain} \
                                            ${residue_entity} \
                                            ${residue_alt} \
                                            ${conf_model} \
	                                    ${small_angle} \
	                                    ${potential} \
					    ${energy_cutoff_atom} \
					    ${energy_cutoff_residue} \
					    ${pdbx_dump_file}
