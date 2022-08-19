#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/double-dimer-library-002.cif

rotag_predict -S ${pdbx_file}
