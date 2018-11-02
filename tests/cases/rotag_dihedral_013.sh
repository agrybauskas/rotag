#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-library-selected-001.cif

rotag_dihedral --radians ${pdbx_file}
