#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/5svd.dump
selection_command="resname ASP & atoname CA"

"$(dirname "$0")"/../scripts/command_line_parser "${selection_command}" \
		                                 "${pdbx_dump_file}"
