#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/aspartic-acid-001.cif

rotag_energy -r ${pdbx_file}
