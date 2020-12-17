#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

atom_names_file=$(dirname "$0")/../inputs/atom-names-003.dat

$(dirname "$0")/../scripts/sort_atom_names ${atom_names_file} 'gn'
