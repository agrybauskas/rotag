#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-001.cif

rotag_scan --parameters 'cutoff_atom=Inf' -H -b '!;OD1.MG=1.050' ${pdbx_file}
