#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

parameters_file=$(dirname "$0")/../inputs/parameters/Parameters-001.cif

$(dirname "$0")/../scripts/obtain_force_field_data ${parameters_file}
