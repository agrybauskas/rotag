#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/phenylalanine-001.cif

rotag_scan ${pdbx_file} 2>&1
