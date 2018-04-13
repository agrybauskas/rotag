#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/5svd.cif
categories="_atom_site,_citation"

$(dirname "$0")/../scripts/obtain_category_data ${categories} ${pdbx_file}
