#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/k-with-sidechains-with-connections-library-002.cif

rotag_add -k -S ${pdbx_file}
