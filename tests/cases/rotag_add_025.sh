#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-library-001.cif

rotag_add -S ${pdbx_file}
