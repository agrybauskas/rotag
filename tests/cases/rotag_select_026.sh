#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/fragment-001.cif

rotag_select --tags '_atom_site' -F 'csv' ${pdbx_file}
