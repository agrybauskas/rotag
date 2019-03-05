#!/bin/bash

atom_coord=$(dirname "$0")/../inputs/atom-coord-001.dat

$(dirname "$0")/../bin/find_euler_angle ${atom_coord}
