#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/polipeptides/dipeptide-001.dump
calc_mainchain=1

$(dirname "$0")/../scripts/all_dihedral ${pdbx_dump_file} $calc_mainchain
