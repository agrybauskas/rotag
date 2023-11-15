#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.cif

rotag_scan -H -a 36.0 ${pdbx_file}
