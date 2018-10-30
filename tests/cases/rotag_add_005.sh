#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/synthetic/xaa-002.cif

rotag_add -H ${pdbx_file}
