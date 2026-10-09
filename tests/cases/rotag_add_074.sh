#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/aspartic-acid-library-003.cif

rotag_add -S ${pdbx_file}
