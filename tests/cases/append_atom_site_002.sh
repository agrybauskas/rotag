#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file_dump_1=$(dirname "$0")/../inputs/amino-acids/valine-001.dump
pdbx_file_dump_2=$(dirname "$0")/../inputs/amino-acids/serine-001.dump

$(dirname "$0")/../scripts/append_atom_site 1 ${pdbx_file_dump_1} ${pdbx_file_dump_2}
