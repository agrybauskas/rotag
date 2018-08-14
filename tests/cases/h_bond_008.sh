#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/synthetic/h-donor-acceptor-008.dump
donor_id=2
acceptor_id=3

$(dirname "$0")/../scripts/h_bond ${donor_id} \
                                  ${acceptor_id} \
                                  ${pdbx_dump_file}
