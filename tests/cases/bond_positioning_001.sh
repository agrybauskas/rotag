#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/alanine-001.dump

$(dirname "$0")/../scripts/bond_positioning ${pdbx_dump_file}
