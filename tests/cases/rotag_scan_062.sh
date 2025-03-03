#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/2qty-multiple-mg-001.cif

rotag_scan -H -b '..120.0..;*.MG=1.2..0.5..3.2;*-*-*.MG=-180.0..18.0..180.0;*-*.MG=81.1..24.7..179.9' ${pdbx_file}
