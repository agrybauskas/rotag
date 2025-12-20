#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-and-mg-with-sidechains-with-connections-002.cif

rotag_angle -H -S ${pdbx_file}
