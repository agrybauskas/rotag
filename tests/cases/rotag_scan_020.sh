#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-library-006.cif

rotag_scan ${pdbx_file} 2>&1
