#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-library-001.cif

rotag_add -r ${pdbx_file}
