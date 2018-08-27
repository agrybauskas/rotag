#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/isoleucine-selected-001.cif

rotag_dihedral ${pdbx_file}
