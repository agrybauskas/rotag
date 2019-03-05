#!/bin/bash

atom_coord=$(dirname "$0")/../inputs/matrices/matrix-001.dat

$(dirname "$0")/../bin/transpose ${atom_coord}
