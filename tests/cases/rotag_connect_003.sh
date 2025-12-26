#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-and-mg-with-sidechains-002.cif

rotag_connect --covale '4,3,6' --hydrog '2,3;4,5' ${pdbx_file}
