#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-001.dump
include_specifier=""
exclude_specifier=""
data_specifier=""
order_specifier=
group_specifier="type_symbol C,N,O"

"$(dirname "$0")"/../scripts/filter "${include_specifier}" \
		                    "${exclude_specifier}" \
                		    "${data_specifier}"    \
		                    "${pdbx_dump_file}"    \
				    "${order_specifier}"   \
				    "${group_specifier}"   \
