#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-011.cif

rotag_scan --angles '-180.0..36.0..180.0' -r ${pdbx_file}
