#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file="$(dirname "$0")/../inputs/synthetic/single-atoms-001.dump"
atom_i_id=1
atom_j_ids="2,3,4"
parameters='lj_k=1.0,c_k=1.0,h_k=1.0,cutoff_atom=0.5,cutoff_residue=1.0,cutoff_start=2.5,cutoff_end=5.0'

$(dirname "$0")/../scripts/lowest_energy_state ${atom_i_id} \
                                               ${atom_j_ids} \
			                                   ${pdbx_dump_file} \
                                               ${parameters}
