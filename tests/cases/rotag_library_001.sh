#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/asparagine-surrounded-001.cif

rotag_library -i ${pdbx_file} --cutoff-atom 20 --cutoff-residue 0
