#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/sidechain-stream-001.cif
item_specifier="_atom_site"

$(dirname "$0")/../scripts/obtain_pdbx_loop ${item_specifier} ${pdbx_file}
