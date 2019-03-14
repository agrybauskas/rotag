#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

parameters_file=$(dirname "$0")/../inputs/parameters/Parameters-001.cif

$(dirname "$0")/../scripts/force_field_parameters ${parameters_file}
