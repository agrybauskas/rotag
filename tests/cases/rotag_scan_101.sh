#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-H-001.cif

rotag_scan --rand-seed 23 --rand-sweep 10 ${pdbx_file}
