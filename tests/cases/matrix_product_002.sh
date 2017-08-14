#!/bin/bash
cd "$(dirname "$0")"

matrix_file=../inputs/matrices_002.dat
variable_values=""

../programs/matrix_product "${variable_values}" ${matrix_file}
