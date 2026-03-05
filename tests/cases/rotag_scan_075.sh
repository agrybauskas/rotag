#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.cif

rotag_scan -H -b '!;OD1.MG=1.050' tests/inputs/hetatoms/mg-with-sidechains-with-connections-001.cif
