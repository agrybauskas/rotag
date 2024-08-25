#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/2bama-active-site-with-ca.cif

rotag_scan -H -a '120.0,OD1.H1=3.0..1.0..4.0' ${pdbx_file}
