#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/threonine-001.cif

rotag_dihedral -S ${pdbx_file}
