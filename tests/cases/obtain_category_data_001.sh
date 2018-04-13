#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/amino-acids/serine-001.cif
categories="_atom_site,_rotamer"

$(dirname "$0")/../scripts/obtain_category_data ${categories} ${pdbx_file}
