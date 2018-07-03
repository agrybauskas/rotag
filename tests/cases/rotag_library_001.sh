#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-001.cif

rotag_library -i ${pdbx_file} --cutoff-atom 10 --cutoff-residue 0
