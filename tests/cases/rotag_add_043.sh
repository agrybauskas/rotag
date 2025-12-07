#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.cif

rotag_add -l 'MG,MG' ${pdbx_file}
