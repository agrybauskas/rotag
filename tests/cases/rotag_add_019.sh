#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-library-007.cif

rotag_add -S ${pdbx_file}
