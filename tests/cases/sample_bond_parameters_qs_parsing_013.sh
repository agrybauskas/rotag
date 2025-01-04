#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

$(dirname "$0")/../scripts/sample_bond_parameters_qs_parsing 'CG-OD1.MG=-90.0..40.0..90.0'
