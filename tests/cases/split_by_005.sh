#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-with-connections-library-001.dump
split_by_attributes="pdbx_PDB_model_num,label_asym_id"
append_to_alt_ids=0
append='label_alt_id=.,pdbx_auth_alt_id=.'

"$(dirname "$0")"/../scripts/split_by ${split_by_attributes} \
                                      ${pdbx_dump_file} \
                                      ${append_to_alt_ids} \
                                      ${append}
