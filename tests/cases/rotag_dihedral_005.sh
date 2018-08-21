#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/phenylalanine-selected-001.cif

rotag_dihedral -r ${pdbx_file}
