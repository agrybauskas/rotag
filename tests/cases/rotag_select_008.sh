#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/lysine-001.cif

rotag_select -t 'atomid 8642 group 3; atomid 8645 group 2' ${pdbx_file}
