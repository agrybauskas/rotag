#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/active-site-003.cif

rotag_scan -H --rand-seed 23 --rand-step 5 --limit 5 ${pdbx_file}
