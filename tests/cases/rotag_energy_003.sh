#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/serine-H-bonding-001.cif

rotag_energy --potential soft_sphere ${pdbx_file}
