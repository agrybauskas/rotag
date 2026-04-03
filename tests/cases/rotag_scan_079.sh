#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-H-001.cif

rotag_scan -x 23 -X 10 ${pdbx_file}
