#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

unique_residue_keys='27,A,1,.,27,A,?;.,B,1,.,101,A,1'

"$(dirname "$0")"/../scripts/group_unique_residue_keys ${unique_residue_keys}
