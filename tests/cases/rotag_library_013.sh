#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-with-connections-library-003.cif

rotag_library -M 0.5 ${pdbx_file}
