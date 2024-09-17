#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/ca-with-sidechains-with-connections-001.cif

rotag_scan -H -a '!; CB-CG-OD1.CA=90.0' ${pdbx_file}
