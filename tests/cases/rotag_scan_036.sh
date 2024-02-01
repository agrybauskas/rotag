#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.cif

rotag_scan -H -a "120.0,CB-CG-OD1.MG=90.0,OD1.MG=0.5..0.5..1.5" ${pdbx_file}
