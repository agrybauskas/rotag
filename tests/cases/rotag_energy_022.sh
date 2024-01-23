#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-library-008.cif

rotag_energy -S -P ${pdbx_file}
