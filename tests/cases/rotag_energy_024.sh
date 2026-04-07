#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/fragments/glutamic-acids-with-magnesium-001.cif

rotag_energy ${pdbx_file}
