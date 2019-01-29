#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

# $(dirname "$0")/../scripts/sample_angles "0..2" 0.2
# $(dirname "$0")/../scripts/sample_angles "0..1" 0.1
$(dirname "$0")/../scripts/sample_angles "-1..1" 0.1
