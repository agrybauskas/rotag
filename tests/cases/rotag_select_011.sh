#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/amino-acids/leucine-selected-001.cif

rotag_select -t 'target && atomname CA' ${pdbx_file}
