#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_004.dat

../programs/matrix_product < ${matrices}
