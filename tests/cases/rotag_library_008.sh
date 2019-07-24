#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/lysine-library-001.cif

rotag_library -t 1 --tags '_[local]_rotamer_angle' -F csv ${pdbx_file}
