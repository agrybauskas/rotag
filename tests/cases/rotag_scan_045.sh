#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-005.cif

rotag_scan -H -a '!; O.MG=3.0..0.5..4.5' ${pdbx_file} 2>&1
