#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/isoleucine-001.cif

cat ${pdbx_file} | rotag_select -t 'atomname CA' -k
