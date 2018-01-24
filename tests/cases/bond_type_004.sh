#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-connected-001.dump

$(dirname "$0")/../scripts/bond_type ${pdbx_dump_file} | sort -k 1 -n
