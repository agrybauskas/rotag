#!/bin/bash
cd "$(dirname "$0")"

two_coord_sets=../inputs/two_coord_sets_001.dat

../programs/rmsd < ${two_coord_sets}
