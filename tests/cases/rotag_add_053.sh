#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-with-connections-library-006.cif

rotag_add -S --optimise-ligands ${pdbx_file}
