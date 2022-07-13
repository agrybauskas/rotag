#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-with-connections-002.dump
atom_id=1925
do_bond_torsion=0
do_bond_stretching=1
do_angle_bending=0

$(dirname "$0")/../scripts/rotation_translation ${atom_id} \
                                                ${pdbx_dump_file} \
                                                ${do_bond_torsion} \
                                                ${do_bond_stretching} \
                                                ${do_angle_bending}
