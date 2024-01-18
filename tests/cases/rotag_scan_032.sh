#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.cif

rotag_scan -H -a '90.0,CB-CG-OD1-MG=90.0' ${pdbx_file}
