#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/arginine-001.cif

rotag_angle -S ${pdbx_file}
