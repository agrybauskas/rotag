#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/cysteine-001.cif

rotag_connect --disulf '1,2' ${pdbx_file}
