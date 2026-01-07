#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/leucine-selected-001.cif

rotag_mutate -m '1:TRP' ${pdbx_file}
