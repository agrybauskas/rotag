#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/sidechain-stream-001.cif

$(dirname "$0")/../scripts/obtain_atom_sites ${pdbx_file}
