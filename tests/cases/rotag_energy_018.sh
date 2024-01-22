#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.cif

rotag_energy -S ${pdbx_file}
