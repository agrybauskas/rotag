#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/5svd.cif
item_specifier="_atom_site,_pdbx_validate_peptide_omega"

$(dirname "$0")/../scripts/obtain_pdbx_loop ${item_specifier} ${pdbx_file}
