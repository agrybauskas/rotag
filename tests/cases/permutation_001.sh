#!/bin/bash
cd "$(dirname "$0")"

set_file=../inputs/set_of_angles_001.dat

../programs/permutation ${set_file}
