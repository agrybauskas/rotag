#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/active-site-library-001.cif

rotag_add -S --dbscan-min 5 --dbscan-distance 2.0 --optimise-ligands ${pdbx_file}
