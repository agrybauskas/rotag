#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/pentamer-library-001.cif

rotag_predict -S ${pdbx_file}
