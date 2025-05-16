#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-library-002.cif

rotag_scan --use-library -H -b '*.MG=2.2..1.0..3.2;*-*-*.MG=-180.0..72.0..180.0;*-*.MG=81.1..49.4..179.9' ${pdbx_file}
