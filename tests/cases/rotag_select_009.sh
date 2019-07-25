#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/models/aspartic-acid-model-001.cif

rotag_select -t 'altid 2' --related-data ${pdbx_file}
