#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/isoleucine-001.cif

rotag_length -S ${pdbx_file}
