#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/5svd.cif
item_specifier="_entry.id,_audit_conform.dict_name,_citation.title,_citation.id,_entity_poly.pdbx_seq_one_letter_code"

$(dirname "$0")/../scripts/obtain_pdbx_line_new ${item_specifier} ${pdbx_file}
