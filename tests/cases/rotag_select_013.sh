#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/3d1g.pdb

rotag_select --pdb -t 'atomid 10' ${pdbx_file}
