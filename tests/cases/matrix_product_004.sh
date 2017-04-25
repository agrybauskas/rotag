#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_004.dat
variable_values="x=36,y=9,z=4"

../programs/matrix_product "${variable_values}" < ${matrices}
