#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-011.cif

rotag_scan -r ${pdbx_file}
