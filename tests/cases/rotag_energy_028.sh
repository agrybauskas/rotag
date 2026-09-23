#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/glutamic-acid-with-mg-003.cif

rotag_energy -S ${pdbx_file}
