#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/proline-002.dump
only_covalent_radii=1

$(dirname "$0")/../scripts/connect_atoms ${pdbx_dump_file} ${only_covalent_radii} \
    | sort -k 1 -n
