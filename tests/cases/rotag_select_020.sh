#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/polipeptides/dipeptide-001.cif

rotag_select -t 'expand (atomname N && chain A && resname CYS)' ${pdbx_file}
