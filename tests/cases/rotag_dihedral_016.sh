#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/glycine-selected-001.cif

rotag_dihedral --radians ${pdbx_file} 2>&1
