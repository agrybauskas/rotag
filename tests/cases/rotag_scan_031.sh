#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/threonine-003.cif

rotag_scan ${pdbx_file} 2>&1
