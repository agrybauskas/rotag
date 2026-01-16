#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-001.cif

rotag_mutate -m '1:CA' ${pdbx_file}
