#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/ca-with-sidechains-001.cif

rotag_mutate -m '2:CA' ${pdbx_file}
