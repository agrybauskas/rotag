#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/serine-001.cif

rotag_scan -M 0.1 --limit 10 ${pdbx_file}
