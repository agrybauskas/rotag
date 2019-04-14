#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/3d1g.pdb

$(dirname "$0")/../scripts/obtain_pdb_atom_site ${pdbx_file}
