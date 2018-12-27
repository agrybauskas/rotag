#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/glutamine-003.cif

rotag_energy -r ${pdbx_file}
