#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_file_dump=$(dirname "$0")/../inputs/amino-acids/arginine-001.dump

$(dirname "$0")/../scripts/prepare_for_calc ${pdbx_file_dump}
