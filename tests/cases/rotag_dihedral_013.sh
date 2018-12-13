#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/polipeptides/dipeptide-001.cif

rotag_dihedral -m ${pdbx_file}
