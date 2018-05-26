#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/lysine-001.dump
conf_model="rotation_only"
potential="leonard_jones"
energy_cutoff_atom=0.04
energy_cutoff_summed="Inf"
residue_id="572"
small_angle="1*pi"

$(dirname "$0")/../scripts/generate_library ${residue_id} \
                              	            ${conf_model} \
	                                    ${small_angle} \
	                                    ${potential} \
					    ${energy_cutoff_atom} \
					    ${energy_cutoff_summed} \
					    ${pdbx_dump_file}
