#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/3qek.cif

rotag_select -t 'id 1' ${pdbx_file}
