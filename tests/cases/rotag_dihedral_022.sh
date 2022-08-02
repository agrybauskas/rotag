#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-with-bb-001.cif

rotag_dihedral -H -S -M ${pdbx_file}
