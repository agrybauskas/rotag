#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-004.cif

rotag_scan -H -a 'N.MG=1.0..0.5..1.5' ${pdbx_file}
