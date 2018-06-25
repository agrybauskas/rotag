#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.dump
add_only_clear_positions=1

$(dirname "$0")/../scripts/add_hydrogens ${pdbx_dump_file} \
                                         ${add_only_clear_positions}
