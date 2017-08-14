#!/bin/bash
cd "$(dirname "$0")"

matrix_file=../inputs/matrices_001.dat
variable_values="x=6,y=3,z=2"

../programs/matrix_product "${variable_values}" ${matrix_file}
