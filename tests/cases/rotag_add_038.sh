#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/na-with-sidechains-with-connections-library-002.cif

rotag_add -S ${pdbx_file}
