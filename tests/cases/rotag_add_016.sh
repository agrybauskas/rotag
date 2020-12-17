#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-library-003.cif

rotag_add -S ${pdbx_file} 2>&1
