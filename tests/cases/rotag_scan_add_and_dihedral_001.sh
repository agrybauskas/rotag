#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.cif

rotag_scan --angles 90.0 ${pdbx_file} | rotag_add -S | rotag_dihedral -S
