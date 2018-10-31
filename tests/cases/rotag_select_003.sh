#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-012.cif

rotag_select -t 'model 2' ${pdbx_file}
