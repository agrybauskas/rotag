#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/zn-surrounded-with-connections-001.cif

rotag_length -H -S -M ${pdbx_file}
