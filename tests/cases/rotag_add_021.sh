#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-with-sidechains-002.cif

rotag_add -H ${pdbx_file}
