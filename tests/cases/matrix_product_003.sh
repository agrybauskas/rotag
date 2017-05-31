#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_003.dat
variable_values=""

../programs/matrix_product "${variable_values}" ${matrices}
