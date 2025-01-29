#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/asn-with-ca-with-connections-001.cif

rotag_scan -H -a '..120.0..' ${pdbx_file}
