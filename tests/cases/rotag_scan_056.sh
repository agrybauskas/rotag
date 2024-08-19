#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.cif

rotag_scan -H -a "CG-OD1.MG=100..36.0..172.0" ${pdbx_file}
