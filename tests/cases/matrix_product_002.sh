#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_005.dat

../programs/matrix_product < ${matrices}
