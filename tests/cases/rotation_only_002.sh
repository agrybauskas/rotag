#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/aspartic_acid_006.cif
target_atom="label_atom_id OD1"

num_of_points=289 # Number of pseudo-atoms that will be generated.

../programs/rotation_only "${target_atom}" "${num_of_points}" < ${cif_file}
