#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_002.dat

../programs/matrix_product < ${matrices}
