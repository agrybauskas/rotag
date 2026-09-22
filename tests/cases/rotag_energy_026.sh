#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-with-mg-library-002.cif

rotag_energy -S ${pdbx_file}
