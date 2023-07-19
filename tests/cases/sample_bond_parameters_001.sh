#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 0
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 1
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 9
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 10
$(dirname "$0")/../scripts/sample_bond_parameters "4..5" 4
$(dirname "$0")/../scripts/sample_bond_parameters "4..5" 5
$(dirname "$0")/../scripts/sample_bond_parameters "0..1" 9
$(dirname "$0")/../scripts/sample_bond_parameters "0..1" 10
$(dirname "$0")/../scripts/sample_bond_parameters "-1..1" 19
$(dirname "$0")/../scripts/sample_bond_parameters "-1..1" 20
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 10 0 0
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 10 0 1
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 10 1 0
$(dirname "$0")/../scripts/sample_bond_parameters "0..2" 10 1 1
