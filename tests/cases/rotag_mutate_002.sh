#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/arginine-selected-002.cif

rotag_mutate -m 'SER' ${pdbx_file}
