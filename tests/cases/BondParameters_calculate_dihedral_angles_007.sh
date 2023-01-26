#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/polipeptides/dipeptide-001.dump
calc_mainchain=1

$(dirname "$0")/../scripts/BondParameters_calculate_dihedral_angles ${pdbx_dump_file} $calc_mainchain
