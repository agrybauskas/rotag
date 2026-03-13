#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/active-site-library-001.cif

rotag_add -S --optimise-ligands ${pdbx_file}
