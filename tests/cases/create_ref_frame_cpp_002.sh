#!/bin/bash

atom_coord=$(dirname "$0")/../inputs/atom-coord-002.dat

$(dirname "$0")/../bin/create_ref_frame ${atom_coord}
