#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file=$(dirname "$0")/../inputs/cysteine-serine-stream-001.cif

$(dirname "$0")/../scripts/split_pdbx_file ${pdbx_file}
