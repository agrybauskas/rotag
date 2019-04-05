#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/arginine-002.dump
potential="composite"

$(dirname "$0")/../scripts/energy ${potential} ${pdbx_dump_file}
