#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/asparagine-001.cif

rotag_select -t 'atomname C' -s 'atomname N' ${pdbx_file}
