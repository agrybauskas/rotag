#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-and-h2o-with-sidechains-with-connections-001.cif

rotag_scan -H -a '120.0,CG-OD1.MG.O=90.0' ${pdbx_file} | rotag_add -S
