#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/lysine-001.cif

rotag_select -t '3:atomid 8642; 2:atomid 8645' ${pdbx_file}
