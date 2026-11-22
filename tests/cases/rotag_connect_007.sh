#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-and-mg-with-sidechains-003.cif

rotag_connect --metalc '1,2' ${pdbx_file}
