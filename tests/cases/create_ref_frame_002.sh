#!/bin/bash
cd "$(dirname "$0")"

atom_coord=../inputs/three_atom_coord_002.dat

../programs/create_ref_frame ${atom_coord}
