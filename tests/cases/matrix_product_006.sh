#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_006.dat
variable_values="x=1,y=1,z=1,chi=0"

../programs/matrix_product "${variable_values}" ${matrices}
