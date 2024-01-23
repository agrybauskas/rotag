#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/serine-with-mg-library-001.cif

rotag_energy -S -P ${pdbx_file}
