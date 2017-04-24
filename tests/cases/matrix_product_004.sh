#!/bin/bash
cd "$(dirname "$0")"

matrices=../inputs/matrices_001.dat

../programs/matrix_product < ${matrices}
