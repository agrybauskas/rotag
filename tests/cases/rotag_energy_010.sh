#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/empty.cif

rotag_energy -S ${pdbx_file} 2>&1
