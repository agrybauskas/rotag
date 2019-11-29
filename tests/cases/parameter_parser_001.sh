#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

parameter_command="resname ASP && atomname CA"

"$(dirname "$0")"/../scripts/parameter_parser "${paramter_command}"
