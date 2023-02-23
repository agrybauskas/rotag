#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/histidine-bond-angle-only-001.dump
residue_unique_key="94,A,1,."
angle_values="CB-CG-ND1 0.6667*pi  & CB-CG-CD2 0.6667*pi & CG-ND1-CE1 0.6111*pi & CG-CD2-NE2 0.6111*pi"
do_angle_bending=1
do_bond_torsion=0
do_bond_stretching=0

$(dirname "$0")/../scripts/generate_rotamer "${residue_unique_key}" \
	                                        "${angle_values}" \
					                        "${pdbx_dump_file}" \
                                            "${do_angle_bending}" \
                                            "${do_bond_torsion}" \
                                            "${do_bond_stretching}"
