#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/methionine-selected-001.cif

rotag_select -t 'target' ${pdbx_file} 2>&1
