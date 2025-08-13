#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/cd-with-sidechains-with-connections-library-001.cif

rotag_add -k -S ${pdbx_file}
