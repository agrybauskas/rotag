#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/h2o-and-mg-with-sidechains-with-connections-002.cif

rotag_scan -H -a '!; OD1.H1-O=..36.0..; H1-O.MG=..36.0..' ${pdbx_file}
