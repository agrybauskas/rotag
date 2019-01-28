#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

$(dirname "$0")/../scripts/sample_angles_qs_parsing 'chi1=0.0..36.0..360.0'
$(dirname "$0")/../scripts/sample_angles_qs_parsing 'chi1=0.0..36.0..360.0,chi2=0.0..36.0..360.0'
$(dirname "$0")/../scripts/sample_angles_qs_parsing '0.0..36.0..360.0'
$(dirname "$0")/../scripts/sample_angles_qs_parsing '0.0..360.0'
