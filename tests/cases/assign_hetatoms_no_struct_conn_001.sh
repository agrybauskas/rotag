#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_struct_conn_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-without-connections-and-struct-conn-001.dump
do_bond_stretching=0
do_angle_bending=1
do_bond_rotation=0

$(dirname "$0")/../scripts/assign_hetatoms_no_struct_conn \
    ${pdbx_struct_conn_dump_file} ${do_bond_stretching} ${do_angle_bending} ${do_bond_rotation} \
    | sort -k 1 -n
