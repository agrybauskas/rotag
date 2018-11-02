#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/arginine-002.cif

rotag_add -H ${pdbx_file}
