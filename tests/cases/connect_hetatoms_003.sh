#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_struct_conn_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-without-connections-and-struct-conn-002.dump

$(dirname "$0")/../scripts/connect_hetatoms \
    ${pdbx_struct_conn_dump_file} \
    | sort -k 1 -n
