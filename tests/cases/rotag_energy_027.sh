#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/glutamic-acid-with-mg-001.cif

rotag_energy -S ${pdbx_file}
