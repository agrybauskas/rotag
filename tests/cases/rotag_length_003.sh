#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-002.cif

rotag_length -H -S -M ${pdbx_file}
