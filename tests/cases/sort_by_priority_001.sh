#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

atom_names_file=$(dirname "$0")/../inputs/atom-names-001.dat

$(dirname "$0")/../scripts/sort_by_priority ${atom_names_file}
