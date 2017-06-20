#!/bin/bash
cd "$(dirname "$0")"

parameter_file=../../parameters/rotatable_bonds.csv

../programs/rotatable_bonds ${parameter_file}
