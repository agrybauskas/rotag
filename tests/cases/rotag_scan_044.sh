#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-003.cif

rotag_scan -H -a '120.0,N.MG=1.0..0.5..1.5,CB-CG-OD1.MG=90.0' ${pdbx_file}
