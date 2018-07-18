#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/fragment-001.cif

rotag_select -t 'byres atomid 2' ${pdbx_file}
