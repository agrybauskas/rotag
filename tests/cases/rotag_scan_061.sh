#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/asn-with-ca-with-connections-001.cif

rotag_scan -H -a '..120.0..;*.CA=1.7..0.4..3.3;*.ZN=1.7..0.3..2.9;*.NA=1.3..0.6..3.7;*-*-*.*=-180.0..18.0..180.0;*-*.*=81.1..24.7..179.9' ${pdbx_file}
