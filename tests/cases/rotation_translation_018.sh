#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/lysine-001.dump
atom_id=8649
calc_hetatoms=0
do_bond_torsion=1
do_bond_stretching=0
do_angle_bending=0

$(dirname "$0")/../scripts/rotation_translation ${atom_id} \
                                                ${pdbx_dump_file} \
                                                ${calc_hetatoms} \
                                                ${do_bond_torsion} \
                                                ${do_bond_stretching} \
                                                ${do_angle_bending}
