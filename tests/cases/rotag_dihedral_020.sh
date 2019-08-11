#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/serine-001.cif

rotag_dihedral ${pdbx_file} 2>&1
