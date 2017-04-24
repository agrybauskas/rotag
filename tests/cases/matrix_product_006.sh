#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_003.dat

../programs/matrix_product < ${matrices}
