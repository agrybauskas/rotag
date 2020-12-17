#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/alanine-001.cif

rotag_library ${pdbx_file} 2>&1
