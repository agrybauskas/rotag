#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/empty.cif

rotag_add -S ${pdbx_file} 2>&1
