#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/proline-001.dump

$(dirname "$0")/../scripts/connection_digest ${pdbx_dump_file}
