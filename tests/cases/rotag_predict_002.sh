#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/asparagine-dimer-library-001.cif

rotag_predict -S ${pdbx_file}
