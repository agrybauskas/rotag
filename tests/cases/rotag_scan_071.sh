#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-and-mg-with-sidechains-with-connections-002.cif

rotag_scan -H -a '!; CB-CG-OD1.H1=170.0..10.0..180.0; CG-OD1.H1-O=170.0..10.0..180.0; OD1.H1-O.MG=170.0..10.0..180.0' ${pdbx_file}
