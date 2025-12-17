#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-and-mg-with-sidechains-001.cif

rotag_connect --hydrog '2,3;4,5' ${pdbx_file}
