#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/arginine-002.dump
residue_unique_key="315,A,1,A,?,?"
moiety="SER"

$(dirname "$0")/../scripts/replace_with_moiety "${residue_unique_key}" \
                                               ${moiety} \
					                           ${pdbx_dump_file}
