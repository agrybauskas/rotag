#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/polipeptides/dipeptide-001.dump
include_mainchain=1
show_values=1

$(dirname "$0")/../scripts/bendable_angles \
               ${pdbx_dump_file} \
               ${include_mainchain} \
               ${show_values}
