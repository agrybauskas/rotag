#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/zn-surrounded-without-connections-001.cif

rotag_angle -H -S -M ${pdbx_file}
