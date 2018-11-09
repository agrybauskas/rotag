#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/polipeptides/dipeptide-001.cif

rotag_add -H ${pdbx_file}
