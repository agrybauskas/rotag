#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/na-with-sidechains-with-connections-library-003.cif

rotag_add -k -S ${pdbx_file}
