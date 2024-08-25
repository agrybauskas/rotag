#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/histidine-bond-length-only-001.dump
residue_unique_key="94,A,1,.,139,A"
angle_values="CA-CB 3.0 & CB-CG 1.0 & CG-ND1 2.5 & CG-CD2 2.0"
do_angle_bending=0
do_bond_torsion=0
do_bond_stretching=1

$(dirname "$0")/../scripts/generate_rotamer "${residue_unique_key}" \
	                                        "${angle_values}" \
					                        "${pdbx_dump_file}" \
                                            "${do_angle_bending}" \
                                            "${do_bond_torsion}" \
                                            "${do_bond_stretching}"