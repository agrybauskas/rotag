#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_struct_conn_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-with-connections-and-struct-conn-001.dump
no_connection_list=0

$(dirname "$0")/../scripts/connect_hetatoms \
    ${pdbx_struct_conn_dump_file} ${no_connection_list} \
    | sort -k 1 -n
