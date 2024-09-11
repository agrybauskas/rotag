#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-005.cif

rotag_scan -H -a 'C-O.MG=..90.0..' ${pdbx_file}
