#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/aspartic-acid-library-001.cif

rotag_add -S -k ${pdbx_file} 2>&1
