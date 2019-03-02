#!/bin/bash

atom_coord=$(dirname "$0")/../inputs/atom-coord-001.dat

$(dirname "$0")/../bin/create_ref_frame ${atom_coord}
