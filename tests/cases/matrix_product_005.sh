#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_005.dat
variable_values="chi=0,x=1,y=2,z=3"

../programs/matrix_product "${variable_values}" < ${matrices}
