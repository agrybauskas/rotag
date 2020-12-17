#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-library-001.cif

rotag_library -F ? ${pdbx_file} 2>&1
