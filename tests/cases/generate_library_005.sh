#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.dump
conf_model="rotation_only"
potential="leonard_jones"
cutoff=0.04
residue_id="219"
small_angle="0.1*pi"

$(dirname "$0")/../scripts/generate_library ${residue_id} \
                              	            ${conf_model} \
	                                    ${small_angle} \
	                                    ${potential} \
					    ${cutoff} \
					    ${pdbx_dump_file}
