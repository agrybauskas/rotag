#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

$(dirname "$0")/../scripts/sample_angles_qs_parsing '36.0; chi1=-180.0..10.0..180.0; PHE,TYR:chi2=-180.0..72.0..180.0'
