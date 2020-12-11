#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/asparagine-dimer-library-001.cif

rotag_predict -S --print-interaction-matrix ${pdbx_file}
