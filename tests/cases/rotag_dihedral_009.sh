#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/threonine-selected-001.cif

rotag_dihedral -r ${pdbx_file}
