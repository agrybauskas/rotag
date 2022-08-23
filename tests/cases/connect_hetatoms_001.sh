#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-001.dump
struct_conn_dump_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-struct-conn-001.dump
no_connection_list=0

$(dirname "$0")/../scripts/connect_hetatoms \
    ${pdbx_dump_file} ${struct_conn_dump_file} ${no_connection_list} \
    | sort -k 1 -n
