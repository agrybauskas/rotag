#!/bin/bash
cd "$(dirname "$0")"

atom_coord=../inputs/three_atom_coord_001.dat

../programs/find_euler_angles < ${atom_coord}
