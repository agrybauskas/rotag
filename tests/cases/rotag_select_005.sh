#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/5svd.cif

rotag_select -x 22 -t '(atomname CA && resname TRP) rand 1' ${pdbx_file}
