#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

atom_coords_file="$(dirname "$0")"/../inputs/atom-coords-001.dat
ref_frame="local"

"$(dirname "$0")"/../scripts/switch_ref_frame ${atom_coords_file} ${ref_frame}
