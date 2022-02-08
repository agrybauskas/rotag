#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/leucine-selected-001.cif

rotag_mutate -r -m '1:SER,chi1=3.14' ${pdbx_file}
