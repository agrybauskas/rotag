#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/ca-with-sidechains-with-connections-library-001.cif

rotag_add -S ${pdbx_file}
