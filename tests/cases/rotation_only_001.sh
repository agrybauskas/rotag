#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/serine_001.cif
target_atom="label_atom_id OG"

num_of_points=17 # Number of pseudo-atoms that will be generated.

../programs/rotation_only "${target_atom}" "${num_of_points}" < ${cif_file}
