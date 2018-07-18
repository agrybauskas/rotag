#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/serine-H-bonding-001.cif

rotag_energy --potential coulomb ${pdbx_file}
