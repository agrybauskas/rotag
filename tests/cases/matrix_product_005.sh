#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_005.dat
variable_values="chi=0"

../programs/matrix_product "${variable_values}" < ${matrices}
