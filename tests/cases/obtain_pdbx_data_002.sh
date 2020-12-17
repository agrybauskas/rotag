#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/5svd.cif
data_identifier="_audit_author,_audit_conform.dict_name,_entity_poly.pdbx_seq_one_letter_code"

cat ${pdbx_file} | $(dirname "$0")/../scripts/obtain_pdbx_data ${data_identifier}
