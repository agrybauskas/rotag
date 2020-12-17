#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/aspartic-acid-001.cif

rotag_add -s -H ${pdbx_file} 2>&1
