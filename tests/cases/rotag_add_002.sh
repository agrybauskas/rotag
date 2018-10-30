#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-library-002.cif

rotag_add -r ${pdbx_file}
