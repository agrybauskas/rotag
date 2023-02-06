#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-001.dump
include_mainchain=0
do_calculation=1

$(dirname "$0")/../scripts/rotatable_bonds \
               ${pdbx_dump_file} \
               ${include_mainchain} \
               ${do_calculation}
