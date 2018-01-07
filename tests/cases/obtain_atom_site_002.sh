#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/5svd.cif

$(dirname "$0")/../scripts/obtain_atom_site ${pdbx_file}
