#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/tryptophan-bond-angle-only-001.dump
residue_unique_key="276,A,1,.,321,A"
angle_values="CG-CD1-NE1 0.5*pi & CD2-CE3-CZ3 0.5*pi"
do_angle_bending=1
do_bond_torsion=0
do_bond_stretching=0

$(dirname "$0")/../scripts/generate_rotamer "${residue_unique_key}" \
	                                        "${angle_values}" \
					                        "${pdbx_dump_file}" \
                                            "${do_angle_bending}" \
                                            "${do_bond_torsion}" \
                                            "${do_bond_stretching}"
