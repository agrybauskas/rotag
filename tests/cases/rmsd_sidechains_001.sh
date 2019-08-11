#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file_1=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-003.dump
pdbx_dump_file_2=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-004.dump

$(dirname "$0")/../scripts/rmsd_sidechains ${pdbx_dump_file_1} ${pdbx_dump_file_2}
