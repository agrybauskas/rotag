#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

atom_names_file=$(dirname "$0")/../inputs/atom-names-001.dat

$(dirname "$0")/../scripts/sort_atom_names ${atom_names_file} 'n' 2>&1 \
    | sed 's/line\s*[0-9]*.$/line <row>./g' \
    | sed 's/0x[0-9a-f]\{12\}/<hex>/g'
