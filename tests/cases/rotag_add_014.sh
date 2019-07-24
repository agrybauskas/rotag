#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/dihedral-angles/serine-dihedral-angles-001.cif

rotag_add -S ${pdbx_file} 2>&1
