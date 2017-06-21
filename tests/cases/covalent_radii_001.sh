#!/bin/bash
cd "$(dirname "$0")"

parameter_file=../../parameters/covalent_radii.csv

../programs/covalent_radii ${parameter_file}
