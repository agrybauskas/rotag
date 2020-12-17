#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-library-002.cif

rotag_add -S ${pdbx_file} 2>&1
