#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-with-connections-001.dump

$(dirname "$0")/../scripts/remove_connections 1925 149,884 ${pdbx_dump_file} | sort -k 1 -n
