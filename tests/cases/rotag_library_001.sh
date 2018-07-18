#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_library --cutoff-atom 10 --cutoff-residue 0 ${pdbx_file}
