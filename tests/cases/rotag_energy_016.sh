#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/glutamic-acid-001.cif

rotag_energy -S -F csv ${pdbx_file}
