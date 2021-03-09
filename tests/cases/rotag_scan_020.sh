#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-library-006.cif

rotag_scan --angles '-180.0..36.0..180.0' ${pdbx_file} 2>&1
