#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-with-sidechains-with-connections-001.cif

rotag_scan -H -a '120.0,CB-CG-OD1.H1=90.0' ${pdbx_file}
