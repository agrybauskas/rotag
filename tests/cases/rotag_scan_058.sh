#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-and-mg-with-sidechains-with-connections-001.cif

rotag_scan -H -a '!; OD1.H1=1.0..1.0..15.0' ${pdbx_file}
