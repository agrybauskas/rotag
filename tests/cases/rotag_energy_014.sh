#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/serine-001.cif

rotag_energy -F csv ${pdbx_file}
