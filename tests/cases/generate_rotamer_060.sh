#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/tryptophan-bond-angle-only-001.dump
residue_unique_key="276,A,1,."
angle_values="chi1 0.5*pi & chi2 0.5*pi"
do_angle_bending=0
do_bond_torsion=1
do_bond_stretching=0

$(dirname "$0")/../scripts/generate_rotamer "${residue_unique_key}" \
	                                        "${angle_values}" \
					                        "${pdbx_dump_file}" \
                                            "${do_angle_bending}" \
                                            "${do_bond_torsion}" \
                                            "${do_bond_stretching}"
