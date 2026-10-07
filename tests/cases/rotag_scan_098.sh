#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/h2o-with-sidechains-with-connections-library-004.cif

rotag_scan -H -b 'chi1=-120.0,chi2=-90.0,CG-OD1.H1=90.0,OD1.H1=50.0,CB-CG-OD1.H1=-180.0,CG-OD1.H1-O=45.0' ${pdbx_file}
