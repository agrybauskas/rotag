#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

unique_residue_keys='71,A,1,.,71,A,.;71,A,1,A,71,A,.;71,A,1,B,71,A,.'

"$(dirname "$0")"/../scripts/group_unique_residue_keys ${unique_residue_keys}
