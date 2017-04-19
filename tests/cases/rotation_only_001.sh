#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/serine_001.cif
target_atom="label_atom_id OG"

../programs/rotation_only "${target_atom}" < ${cif_file}
