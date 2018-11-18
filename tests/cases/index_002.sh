#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/5svd.dump

"$(dirname "$0")"/../scripts/index ${pdbx_dump_file}
