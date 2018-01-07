#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

matrices_file=$(dirname "$0")/../inputs/matrices/matrices-001.dat

$(dirname "$0")/../scripts/flatten ${matrices_file}
