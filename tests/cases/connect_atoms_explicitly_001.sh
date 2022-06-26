#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/mg-with-sidechains-001.dump

$(dirname "$0")/../scripts/connect_atoms_explicitly 1925 149,884 ${pdbx_dump_file} | sort -k 1 -n
