#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-with-connections-library-001.cif

rotag_energy -S -P ${pdbx_file}
