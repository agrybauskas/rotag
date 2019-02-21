#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/empty.cif

rotag_dihedral -S -M ${pdbx_file} 2>&1
