#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

$(dirname "$0")/../scripts/sample_angles_qs_parsing 'TYR:chi1=-180.0..36.0..180.0,chi2=-180.0..36.0..180.0'
