#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/arginine-library-001.cif

rotag_add -r ${pdbx_file}
