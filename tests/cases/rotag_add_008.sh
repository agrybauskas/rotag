#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetero-atoms-001.cif

rotag_add -H ${pdbx_file} 2>&1
