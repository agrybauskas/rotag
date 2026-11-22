#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/serine-selected-002.cif

rotag_length -F csv ${pdbx_file}
