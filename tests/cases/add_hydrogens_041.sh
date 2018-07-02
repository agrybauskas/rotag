#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/surrounded/dipeptide-001.dump

$(dirname "$0")/../scripts/add_hydrogens ${pdbx_dump_file}
