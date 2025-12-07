#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-003.cif

rotag_connect --metalc '2,3' ${pdbx_file}
