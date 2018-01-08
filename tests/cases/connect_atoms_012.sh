#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/lysine-001.dump

$(dirname "$0")/../scripts/connect_atoms ${pdbx_dump_file} | sort -k 1 -n
