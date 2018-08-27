#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif
pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-library-001.cif

rotag_library -c 'Inf' -C 'Inf' -t 1 ${pdbx_file}
rotag_dihedral -r ${pdbx_file}
rotag_energy -r ${pdbx_file}
