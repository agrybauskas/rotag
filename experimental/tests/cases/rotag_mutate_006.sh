#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/leucine-selected-001.cif

cat ${pdbx_file} | rotag_mutate -m 'SER,chi1=90.0' 2>&1
