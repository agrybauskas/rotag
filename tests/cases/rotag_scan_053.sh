#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-007.cif

rotag_scan -H -a 'C-CA-N.MG=90.0,CA-C-O.MG=90.0' ${pdbx_file}
