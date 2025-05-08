#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-with-connections-007.cif

rotag_scan -H -a '!; N.MG=3.0..0.5..3.5;OD1.MG=3.0..0.5..3.5;O.MG=3.0..0.5..3.5' ${pdbx_file}
