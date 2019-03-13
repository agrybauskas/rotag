#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/amino-acids/serine-001.cif

$(dirname "$0")/../scripts/obtain_atom_site_new ${pdbx_file}
