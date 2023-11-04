#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/trimer-library-002.cif

rotag_predict -S ${pdbx_file}
