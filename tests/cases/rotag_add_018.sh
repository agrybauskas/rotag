#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/polipeptides/dipeptide-001.cif
rotamer_angle_file=$(dirname "$0")/../inputs/polipeptides/dipeptide-001.cif

rotag_add -S -a ${rotamer_angle_file} ${pdbx_file}
