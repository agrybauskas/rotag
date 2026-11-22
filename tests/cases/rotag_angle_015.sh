#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/serine-selected-002.cif

rotag_angle -r ${pdbx_file}
