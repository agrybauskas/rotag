#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/5svd.cif
item_specifier="_entry.id,_citation.pdbx_database_id_DOI"

$(dirname "$0")/../scripts/obtain_pdbx_line ${item_specifier} ${pdbx_file}
