#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-008.dump
residue_unique_key="18,A,1,.,?,?"
angle_values="CA-CB-OG 0.5*pi & chi1 0.25*pi"
do_angle_bending=1
do_bond_torsion=1
do_bond_stretching=1

$(dirname "$0")/../scripts/generate_rotamer "${residue_unique_key}" \
	                                        "${angle_values}" \
					                        "${pdbx_dump_file}" \
                                            "${do_angle_bending}" \
                                            "${do_bond_torsion}" \
                                            "${do_bond_stretching}"
