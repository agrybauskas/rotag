#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/fragments/glutamine-003.dump
split_by_attributes="pdbx_PDB_model_num,label_asym_id,label_alt_id"
append_to_alt_ids=1

"$(dirname "$0")"/../scripts/split_by ${split_by_attributes} \
                                      ${pdbx_dump_file} \
                                      ${append_to_alt_ids}
