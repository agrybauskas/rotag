#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.cif

rotag_scan -H -a "!; OD1.MG=0.5..0.5..1.5" ${pdbx_file}
