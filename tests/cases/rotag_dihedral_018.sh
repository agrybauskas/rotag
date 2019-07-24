#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/leucine-001.cif

cat ${pdbx_file} | rotag_dihedral -S -F '?' 2>&1
