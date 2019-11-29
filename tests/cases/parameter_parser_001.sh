#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

parameter_cmd="b={a=5,c=7,d=8}"

"$(dirname "$0")"/../scripts/parameter_parser ${parameter_cmd}
