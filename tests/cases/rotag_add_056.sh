#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/glutamic-acid-with-mg-library-001.cif

rotag_add -S ${pdbx_file}
