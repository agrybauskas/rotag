#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 0.2
$(dirname "$0")/../scripts/sample_bond_parameters "4..5" 0.2
$(dirname "$0")/../scripts/sample_bond_parameters "0..1" 0.1
$(dirname "$0")/../scripts/sample_bond_parameters "-1..1" 0.1
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 0.1 0 1
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 0.1 1 0
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 0.1 1 1
