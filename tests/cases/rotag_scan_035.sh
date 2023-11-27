#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.cif

rotag_scan -H -a "72.0,OD1-MG=0.5..0.5..1.5,CG-OD1-MG=72.0" ${pdbx_file}
