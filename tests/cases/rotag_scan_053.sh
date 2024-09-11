#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-007.cif

rotag_scan -H -a '..120.0..;C-CA-N.MG=..90.0..;CB-CG-OD1.MG=..90.0..;CA-C-O.MG=..90.0..' ${pdbx_file}
