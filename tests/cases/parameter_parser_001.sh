#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

parameter_cmd="a=[{a=5,b=6},{c=7},{d=8}],b=5"

"$(dirname "$0")"/../scripts/parameter_parser ${parameter_cmd}
