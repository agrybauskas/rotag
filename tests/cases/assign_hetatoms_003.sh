#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_struct_conn_dump_file=$(dirname "$0")/../inputs/amino-acids/h2o-and-mg-with-sidechains-with-connections-and-struct-conn-001.dump

$(dirname "$0")/../scripts/assign_hetatoms \
    ${pdbx_struct_conn_dump_file} \
    | sort -k 1 -n
