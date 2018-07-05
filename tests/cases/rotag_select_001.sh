#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/fragment-001.cif

rotag_select -i ${pdbx_file} -t 'byres atomid 2'
