#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-and-mg-with-sidechains-with-connections-002.cif

rotag_scan -H -a '!; OD1.MG=3.0..0.5..3.5' ${pdbx_file}
