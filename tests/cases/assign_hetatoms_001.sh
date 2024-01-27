#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_struct_conn_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-without-connections-and-struct-conn-002.dump
assume_connections=1

$(dirname "$0")/../scripts/assign_hetatoms \
    ${pdbx_struct_conn_dump_file} ${assume_connections} \
    | sort -k 1 -n
