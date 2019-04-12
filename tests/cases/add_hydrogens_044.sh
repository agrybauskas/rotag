#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/fragments/sp3-003.dump

$(dirname "$0")/../scripts/add_hydrogens ${pdbx_dump_file}
