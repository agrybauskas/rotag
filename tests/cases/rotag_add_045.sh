#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/h2o-and-mg-with-sidechains-with-connections-library-003.cif

rotag_add -S ${pdbx_file}
