#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.cif
rotamer_angle_file=$(dirname "$0")/../inputs/libraries/aspartic-acid-library-001.cif

rotag_add -S -a ${rotamer_angle_file} ${pdbx_file}
