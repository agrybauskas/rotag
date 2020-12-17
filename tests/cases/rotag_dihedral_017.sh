#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/leucine-001.cif

rotag_dihedral -S -F csv ${pdbx_file}
