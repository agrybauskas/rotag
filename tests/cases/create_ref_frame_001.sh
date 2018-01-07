#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

atom_coord=$(dirname "$0")/../inputs/atom-coord-001.dat

$(dirname "$0")/../scripts/create_ref_frame ${atom_coord}
