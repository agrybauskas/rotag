#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/arginine-002.dump

$(dirname "$0")/../scripts/energy ${pdbx_dump_file}
