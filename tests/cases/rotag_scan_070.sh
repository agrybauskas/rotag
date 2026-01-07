#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-and-mg-with-sidechains-with-connections-002.cif

rotag_scan -H -a '!; OD1.H1-O=120.0..36.0..156.0; H1-O.MG=120.0..36.0..156.0' ${pdbx_file}
