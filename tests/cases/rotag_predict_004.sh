#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/1l2y.cif

rotag_predict -S --print-adjacency-matrix ${pdbx_file}
