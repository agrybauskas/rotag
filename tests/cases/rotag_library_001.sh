#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/asparagine-surrounded-001.cif

rotag_library -i ${pdbx_file} --cutoff-atom 10 --cutoff-residue 0 --top-rank 1
