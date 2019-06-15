#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/libraries/aspartic-acid-library-001.cif

$(dirname "$0")/../scripts/related_category_data ${pdbx_file}
