#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/2qty-multiple-mg-002.cif

rotag_scan -H -p 'cutoff_atom=115.0' -b '..120.0..;*.MG=2.2..1.0..3.2;*-*-*.MG=-180.0..144.0..180.0;*-*.MG=81.1..49.4..179.9' ${pdbx_file}
