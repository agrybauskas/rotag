#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-H-001.cif

rotag_scan -X 20 ${pdbx_file}
