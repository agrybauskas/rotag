#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/5svd.dump
selection_command="resname ASP && atomname CA"

"$(dirname "$0")"/../scripts/selection_parser "${selection_command}" \
		                              "${pdbx_dump_file}"
