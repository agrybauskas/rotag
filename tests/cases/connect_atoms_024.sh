#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-001.dump
no_connection_list=0
assign_hetatoms=1

$(dirname "$0")/../scripts/connect_atoms \
    ${pdbx_dump_file} ${no_connection_list} ${assign_hetatoms} \
    | sort -k 1 -n
