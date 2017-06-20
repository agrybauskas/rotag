#!/bin/bash
cd "$(dirname "$0")"

parameter_file=../../parameters/vdw_radii.csv

../programs/vdw_radii ${parameter_file}
